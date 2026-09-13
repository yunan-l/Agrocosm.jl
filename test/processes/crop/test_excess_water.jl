using Test
using Agrocosm

isdefined(@__MODULE__, :test_model_state) ||
    include(joinpath(@__DIR__, "..", "..", "helpers", "model_state_fixture.jl"))

# Excess-water damage on an absolute daily rainfall threshold, acting on the
# fraction of the standing crop recovered at harvest.
#
# The load-bearing test here is not the kernel - it is CONSERVATION. This is the
# first mechanism in the project that changes what happens to carbon AT harvest,
# and grain that is not recovered has to stay in the field rather than vanish.

const TW = Float32

function wet_cft(; rate, threshold = 20.0, tolerance = 0.0)
    Agrocosm.CFTParameters{TW, Int32}(;
        (f => (f === :heavy_rain_rate ? TW(rate) :
               f === :heavy_rain_threshold ? TW(threshold) :
               f === :heavy_rain_tolerance ? TW(tolerance) :
               getfield(Agrocosm.cft3, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)
end

@testset "the daily excess is millimetres above the threshold" begin
    ex(mm, thr, growing) = Agrocosm.heavy_rain_excess_today(TW(mm), TW(thr), growing)
    @test ex(30.0, 20.0, true) ≈ TW(10.0)
    @test ex(40.0, 20.0, true) ≈ TW(20.0)      # linear, not a count
    @test ex(20.0, 20.0, true) == zero(TW)     # at the threshold, nothing
    @test ex(5.0, 20.0, true) == zero(TW)      # below it, nothing
    @test ex(60.0, 20.0, false) == zero(TW)    # not growing, nothing
    # A 60 mm day is not a 20.1 mm day, which is the whole reason it is not a
    # day count - and a count is not differentiable either.
    @test ex(60.0, 20.0, true) > ex(20.1, 20.0, true)
end

@testset "the tolerance makes this a wet-YEAR term, not a wet-climate tax" begin
    rec(excess, tol, rate) =
        Agrocosm.excess_water_recovery(TW(excess), TW(tol), TW(rate))
    # A season inside the tolerance loses nothing, however wet the climate.
    @test rec(100.0, 150.0, 0.002) == one(TW)
    @test rec(150.0, 150.0, 0.002) == one(TW)
    # Beyond it, linear in the overshoot.
    @test rec(250.0, 150.0, 0.002) ≈ one(TW) - TW(0.2)
    @test rec(350.0, 150.0, 0.002) ≈ one(TW) - TW(0.4)
    # Clamped, so an extreme season cannot produce a negative harvest.
    @test rec(10_000.0, 150.0, 0.002) == zero(TW)
    # Rate zero is inert whatever the season did - the ablation contract.
    for excess in TW[0.0, 150.0, 900.0]
        @test rec(excess, 150.0, 0.0) == one(TW)
    end
    # This is the failure the first version had: at tolerance zero a median
    # season is charged, which took 54% of rice yield in an AVERAGE year.
    @test rec(116.3, 0.0, 0.002) < rec(116.3, 139.1, 0.002)
    @test rec(116.3, 139.1, 0.002) == one(TW)
end

@testset "the kernel accumulates the season's excess" begin
    crop = test_model_state(init_crop(4, identity))
    state = Agrocosm.crop_prognostic(crop)
    state.phenology.heavy_rain_excess .= zero(TW)
    state.phenology.is_growing .= Int32(1)
    Agrocosm.excess_water!(wet_cft(rate = 0.002), crop, fill(TW(80.0), 4))
    @test all(state.phenology.heavy_rain_excess .≈ TW(60.0))
    Agrocosm.excess_water!(wet_cft(rate = 0.002), crop, fill(TW(80.0), 4))
    @test all(state.phenology.heavy_rain_excess .≈ TW(120.0))
    # A dry day adds nothing, and a stand that is not growing is untouched.
    Agrocosm.excess_water!(wet_cft(rate = 0.002), crop, fill(TW(5.0), 4))
    @test all(state.phenology.heavy_rain_excess .≈ TW(120.0))
    state.phenology.is_growing .= Int32(0)
    Agrocosm.excess_water!(wet_cft(rate = 0.002), crop, fill(TW(200.0), 4))
    @test all(state.phenology.heavy_rain_excess .≈ TW(120.0))
end

@testset "the accumulation is irreversible within a season" begin
    crop = test_model_state(init_crop(1, identity))
    state = Agrocosm.crop_prognostic(crop)
    state.phenology.heavy_rain_excess .= zero(TW)
    state.phenology.is_growing .= Int32(1)
    cft = wet_cft(rate = 0.002)
    for _ in 1:5
        Agrocosm.excess_water!(cft, crop, fill(TW(40.0), 1))
    end
    # A dry spell afterwards does not put the crop back on its feet.
    before = state.phenology.heavy_rain_excess[1]
    for _ in 1:20
        Agrocosm.excess_water!(cft, crop, fill(TW(0.0), 1))
    end
    @test state.phenology.heavy_rain_excess[1] == before
end

@testset "sowing restores the recovery fraction" begin
    # Without the reset, one wet season would follow the field forever - the same
    # reason `grain_set_fraction` is reset in `cultivate!`.
    source = read(joinpath(@__DIR__, "..", "..", "..", "src", "processes", "crop",
                           "cultivate.jl"), String)
    @test occursin("heavy_rain_excess[cell] = zero(T)", source)
    @test :heavy_rain_excess in
          fieldnames(typeof(Agrocosm.crop_prognostic(
              test_model_state(init_crop(1, identity))).phenology))
end

@testset "unrecovered grain stays in the field" begin
    # Conservation. Grain that is not brought in is not destroyed: it must
    # appear as surface litter, and the harvest export must drop by the same
    # amount, or `check_conservation_gates.py` stops closing.
    source = read(joinpath(@__DIR__, "..", "..", "..", "src", "processes", "crop",
                           "harvesting.jl"), String)
    @test occursin("unrecovered_carbon = storage_carbon[cell] * (one(T) - recovery)", source)
    @test occursin("carbon_input[SURFACE_LITTER, cell] = carbon_residue + unrecovered_carbon", source)
    @test occursin("nitrogen_input[SURFACE_LITTER, cell] = nitrogen_residue + unrecovered_nitrogen", source)
    # The export takes the RECOVERED yield, so total out is unchanged.
    @test occursin("carbon_harvest_export[cell] = crop_yield[cell] +", source)
    @test occursin("harvest_nitrogen[cell] = storage_nitrogen[cell] * recovery +", source)
end

@testset "the thresholds are the ones the sweep returned" begin
    @test Agrocosm.cft1.heavy_rain_threshold == TW(10.0)   # wheat, +1.492
    @test Agrocosm.cft2.heavy_rain_threshold == TW(20.0)   # rice, +1.344
    @test Agrocosm.cft3.heavy_rain_threshold == TW(10.0)   # maize, +1.695
    @test Agrocosm.cft9.heavy_rain_threshold == TW(20.0)   # soybean, +1.508
    # Tolerances are the AREA-WEIGHTED 90th percentile of the accumulation each
    # crop meets - weighted because a percentile over cells counts Californian
    # and Mekong rice alike, and the crop is not distributed that way:
    # wet cell-years are 12-14% of disasters, so the damage occupies that decile.
    @test Agrocosm.cft1.heavy_rain_tolerance == TW(179.9)
    @test Agrocosm.cft2.heavy_rain_tolerance == TW(445.2)
    @test Agrocosm.cft3.heavy_rain_tolerance == TW(432.8)
    @test Agrocosm.cft9.heavy_rain_tolerance == TW(271.9)
    for cft in (Agrocosm.cft1, Agrocosm.cft2, Agrocosm.cft3, Agrocosm.cft9)
        # The tolerance must exceed the MEDIAN season or the term is a tax on a
        # wet climate rather than a wet year.
        @test cft.heavy_rain_tolerance > TW(150.0)
    end
    # The rate is now PER CROP, and wheat's is zero. Scoring the recovery
    # fraction against detrended GDHY anomalies per cell, area-weighted, the
    # shipped tolerance gives rho = +0.069 rice, +0.058 maize, +0.008 soybean
    # and -0.042 WHEAT: the wrong sign. Asserted here because a future edit that
    # restores a shared rate would silently re-enable a term measured to be
    # backwards over 162 Mha.
    @test Agrocosm.cft1.heavy_rain_rate == TW(0.0)
    for cft in (Agrocosm.cft2, Agrocosm.cft3, Agrocosm.cft9)
        @test cft.heavy_rain_rate == TW(0.002)
    end
end

@testset "a zero per-crop rate is inert, whatever the season accumulated" begin
    # The ablation contract applied per CROP rather than per arm: with wheat's
    # rate at zero, an `excess_water` run must leave wheat's yield exactly where
    # a run without the mechanism left it, however wet the season was.
    for excess in (TW(0.0), TW(500.0), TW(5000.0))
        @test Agrocosm.excess_water_recovery(
            excess, Agrocosm.cft1.heavy_rain_tolerance,
            Agrocosm.cft1.heavy_rain_rate) == one(TW)
    end
    # And the same season through rice's rate must NOT be inert, or the test
    # above passes because the helper is broken rather than because wheat is off.
    @test Agrocosm.excess_water_recovery(
        TW(5000.0), Agrocosm.cft2.heavy_rain_tolerance,
        Agrocosm.cft2.heavy_rain_rate) < one(TW)
end

@testset "the flag reaches the model through the public entry" begin
    @test :excess_water in Base.kwarg_decl(first(methods(Agrocosm.initialize_simulation)))
    @test :excess_water in Base.kwarg_decl(first(methods(Agrocosm._daily_crop!)))
    @test :excess_water in fieldnames(Agrocosm.SimulationConfiguration)
    configuration = Agrocosm.ablation_excess_water_configuration()
    @test configuration.excess_water === true
    accepted = Set(Base.kwarg_decl(first(methods(Agrocosm.initialize_simulation))))
    for key in keys(configuration)
        @test key in accepted
    end
    @test_throws ArgumentError Agrocosm.ablation_excess_water_configuration(excess_water = false)
end
