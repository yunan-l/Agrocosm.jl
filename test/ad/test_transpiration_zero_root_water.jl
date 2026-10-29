using Agrocosm, Enzyme, Test

function root_water_transpiration_fixture(::Type{T}) where {T <: AbstractFloat}
    return (
        conductance=T[0], assimilation=T[0], co2=T[400], daylength=T[8],
        fpar=T[0], transpiration=zeros(T, 5, 1), demand_sum=T[10],
        supply_sum=T[8], deficit=T[0], sufficiency=T[1], root_carbon=T[20],
        canopy_wet=T[0], is_growing=Int32[1], eeq=T[2.3],
        roots=T[0.5, 0.3, 0.15, 0.04, 0.01], root_water=T[0],
        relative_water=zeros(T, 5, 1), holding_storage=fill(T(50), 5, 1),
    )
end

function root_water_transpiration!(state, cft, parameters)
    Agrocosm.launch_1D!(Agrocosm.water_demand_supply_kernel!,
        state.conductance, state.assimilation, state.co2, state.daylength,
        state.fpar, state.transpiration, state.demand_sum, state.supply_sum,
        state.deficit, state.sufficiency, state.root_carbon, state.canopy_wet,
        state.is_growing, state.eeq, state.roots, state.root_water,
        state.relative_water, state.holding_storage, cft,
        (; lpjmlparams=parameters.lpjml, soil_layers=5, use_precomputed_conductance=true))
    return sum(state.transpiration)
end

for T in (Float32, Float64), crop in (Agrocosm.cft1, Agrocosm.cft3)
    cft = Agrocosm.convert_precision(T, crop)
    parameters = Agrocosm.ModelParameters(T)
    @testset "$T $(cft.path == 1 ? "C3" : "C4") zero root-zone water" begin
        state = root_water_transpiration_fixture(T)
        shadow = Enzyme.make_zero(state)
        result = Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            root_water_transpiration!, Enzyme.Duplicated(state, shadow),
            Enzyme.Const(cft), Enzyme.Const(parameters))
        @test iszero(result[2])
        @test all(iszero, state.transpiration)
        @test all(iszero, state.relative_water)
        @test all(array -> !(eltype(array) <: AbstractFloat) || all(isfinite, array), shadow)

        positive = root_water_transpiration_fixture(T)
        positive.relative_water .= T(0.3)
        positive.conductance .= T(10)
        positive.fpar .= T(0.5)
        positive.eeq .= T(10)
        positive_shadow = Enzyme.make_zero(positive)
        differentiated = deepcopy(positive)
        result = Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.ReverseWithPrimal),
            root_water_transpiration!, Enzyme.Duplicated(differentiated, positive_shadow),
            Enzyme.Const(cft), Enzyme.Const(parameters))
        @test result[2] > zero(T)
        @test all(array -> !(eltype(array) <: AbstractFloat) || all(isfinite, array), positive_shadow)
        projection = sum(positive_shadow.relative_water)
        @test projection > zero(T)
        for step in (T === Float32 ? (T(1e-3), T(5e-4)) : (T(1e-5), T(1e-6)))
            plus = deepcopy(positive)
            minus = deepcopy(positive)
            plus.relative_water .+= step
            minus.relative_water .-= step
            difference = (root_water_transpiration!(plus, cft, parameters) -
                root_water_transpiration!(minus, cft, parameters)) / (2step)
            @test projection ≈ difference rtol=(T === Float32 ? 0.005 : 1e-7)
        end
    end
end
