# A senescent crop can hold exactly zero leaf carbon: LPJmL's own senescence
# clamp sets it, and `senescent_leaf_release` makes it routine. That state used
# to produce a NaN reverse-mode gradient. The cause was
# `nitrogen_demand.jl`'s leaf N:C ratio, whose `leafc > 0 ? n / leafc : 0`
# guard LLVM speculates into a `select` -- see `guarded_quotient`.
#
# Seeds are all zero here, which is what makes this test sharp: with a zero
# cotangent every *finite* local derivative must yield exactly zero, so a
# non-finite entry anywhere in the shadow can only come from a singular local
# derivative. Nothing else in the model has to be involved for this to fail.
@testset "Zero leaf carbon keeps the reverse pass finite" begin
    for T in (Float32, Float64)
        # Precision has to match end to end: a Float32 `cft3` against a Float64
        # state makes the kernel's parameter bundle a mixed-precision union,
        # which Enzyme's type analysis rejects outright.
        cft = Agrocosm.convert_precision(T, cft3)
        parameters = Agrocosm.LPJmLParams{T}()
        crop = init_crop(T, 3, identity)
        state = test_model_state(crop)
        crop.state.phenology.is_growing .= Int32(1)
        crop.state.phenology.senescence .= true
        crop.auxiliary.photosynthesis.lambda .= T(0.4)
        crop.auxiliary.photosynthesis.potential_vcmax .= T(0.00777)
        # The state the release drives a maize stand into: no leaf carbon, but
        # leaf nitrogen and the other organs still stocked.
        crop.state.carbon.leaf .= zero(T)
        crop.state.nitrogen.leaf .= T(0.105)
        crop.state.carbon.root .= T(14.74)
        crop.state.nitrogen.root .= T(0.376)
        crop.state.carbon.storage .= T(48.3)
        crop.state.carbon.pool .= T(55.3)
        crop.state.nitrogen.total .= T(2.49)
        temperature = T[25]

        demand!(st, cft, temp) = begin
            Agrocosm.ndemand_crop!(st, cft,
                Agrocosm.crop_photosynthesis_auxiliary(st).potential_vcmax, temp;
                include_storage_reserve = true, require_active_photosynthesis = true,
                lpjmlparams = parameters)
            nothing
        end

        shadow = Agrocosm.enzyme_zero_tangent(state)
        dtemperature = zeros(T, 1)
        Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), demand!,
            Enzyme.Duplicated(state, shadow), Enzyme.Const(cft),
            Enzyme.Duplicated(temperature, dtemperature))

        @test all(isfinite, dtemperature)
        @test all(isfinite, Agrocosm.crop_prognostic(shadow).carbon.leaf)
        @test all(isfinite, Agrocosm.crop_photosynthesis_auxiliary(shadow).potential_vcmax)
        @test all(isfinite, Agrocosm.crop_stress_auxiliary(shadow).nitrogen_demand_leaf)

        # A tiny positive leaf carbon took the other branch and was always
        # finite; pin it so a future change cannot "fix" the zero case by
        # breaking this one.
        crop.state.carbon.leaf .= T(1e-12)
        tiny_shadow = Agrocosm.enzyme_zero_tangent(state)
        fill!(dtemperature, zero(T))
        Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), demand!,
            Enzyme.Duplicated(state, tiny_shadow), Enzyme.Const(cft),
            Enzyme.Duplicated(temperature, dtemperature))
        @test all(isfinite, dtemperature)
        @test all(isfinite, Agrocosm.crop_prognostic(tiny_shadow).carbon.leaf)
    end
end
