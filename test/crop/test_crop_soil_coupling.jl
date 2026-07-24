using Agrocosm
using Terrarium
using Test

# Conservation of the crop↔soil flux coupling. The crop produces 0D (per-area) litterfall and uptake
# fluxes that the soil biogeochemistry distributes over the root zone as per-volume rates
# (flux·root_fraction/Δz). Because the root fraction sums to unity over the column, the distributed
# input integrated over the column (Σ ·Δz) must recover the original 0D flux — no mass created or
# lost. Each coupling term is isolated with the difference between the tendencies with the crop
# fluxes off and on, so the (identical) decomposition/nitrification baseline cancels exactly.
@testset "Crop → soil flux coupling" begin
    Δz = 0.1
    N = 10
    grid = ColumnGrid(CPU(), UniformSpacing(Δz = Δz, N = N))
    model = CropModel(grid, crop_pft("maize"); soil_hydrology = SoilHydrology(eltype(grid)))
    integrator = initialize(model; initializers = (temperature = 15.0,))
    state = integrator.state
    bgc = model.soil.biogeochem

    # Deterministic soil environment for the decomposition/nitrification baseline.
    set!(state.temperature, 15.0)
    set!(state.saturation_water_ice, 0.5)
    Terrarium.compute_auxiliary!(state, grid, bgc)

    dtend = Terrarium.tendency_fields(state, bgc)
    column(field) = copy(interior(field)[1, 1, :])
    integral(profile) = sum(profile) * Δz   # ∫ · dz over the column, per unit area
    # The root fraction (a normalized profile summing to unity — see test_root_distribution.jl) is what
    # makes the distribution below mass-conserving: the column integral of each input recovers its flux.

    # Baseline tendencies with every crop flux off.
    set!(state.crop_litterfall_carbon, 0.0)
    set!(state.crop_litterfall_nitrogen, 0.0)
    set!(state.crop_nitrogen_uptake, 0.0)
    Terrarium.compute_tendencies!(state, grid, bgc)
    base_litter = column(dtend.litter_carbon)
    base_ammonium = column(dtend.soil_ammonium)
    base_nitrate = column(dtend.soil_nitrate)

    @testset "litterfall carbon enters the litter pool (mass-conserving)" begin
        litterfall_carbon = 1.0e-6   # kgC/m²/s
        set!(state.crop_litterfall_carbon, litterfall_carbon)
        Terrarium.compute_tendencies!(state, grid, bgc)
        Δlitter = column(dtend.litter_carbon) .- base_litter
        @test integral(Δlitter) ≈ litterfall_carbon rtol = 1e-10   # column integral recovers the flux
        @test all(≥(0.0), Δlitter)                                 # every rooted layer receives litter
        set!(state.crop_litterfall_carbon, 0.0)
    end

    @testset "litterfall nitrogen mineralizes to ammonium (mass-conserving)" begin
        litterfall_nitrogen = 3.0e-8   # kgN/m²/s
        set!(state.crop_litterfall_nitrogen, litterfall_nitrogen)
        Terrarium.compute_tendencies!(state, grid, bgc)
        Δammonium = column(dtend.soil_ammonium) .- base_ammonium
        Δnitrate = column(dtend.soil_nitrate) .- base_nitrate
        @test integral(Δammonium) ≈ litterfall_nitrogen rtol = 1e-10   # all of it to ammonium
        @test all(≈(0.0), Δnitrate)                                    # nitrate is untouched
        set!(state.crop_litterfall_nitrogen, 0.0)
    end

    @testset "crop uptake draws down mineral N, split by pool share" begin
        # Equal ammonium and nitrate pools → the uptake splits evenly between them.
        set!(state.soil_ammonium, 0.05)
        set!(state.soil_nitrate, 0.05)
        Terrarium.compute_tendencies!(state, grid, bgc)
        b_ammonium = column(dtend.soil_ammonium)
        b_nitrate = column(dtend.soil_nitrate)

        uptake = 2.0e-8   # kgN/m²/s
        set!(state.crop_nitrogen_uptake, uptake)
        Terrarium.compute_tendencies!(state, grid, bgc)
        Δammonium = column(dtend.soil_ammonium) .- b_ammonium
        Δnitrate = column(dtend.soil_nitrate) .- b_nitrate
        @test integral(Δammonium) ≈ -uptake / 2 rtol = 1e-8   # half drawn from ammonium
        @test integral(Δnitrate) ≈ -uptake / 2 rtol = 1e-8    # half drawn from nitrate
        @test integral(Δammonium .+ Δnitrate) ≈ -uptake rtol = 1e-8   # total uptake conserved
        set!(state.crop_nitrogen_uptake, 0.0)
    end
end
