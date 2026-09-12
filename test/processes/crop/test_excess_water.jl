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

function wet_cft(; rate, threshold = 20.0)
    Agrocosm.CFTParameters{TW, Int32}(;
        (f => (f === :heavy_rain_rate ? TW(rate) :
               f === :heavy_rain_threshold ? TW(threshold) :
               getfield(Agrocosm.cft3, f))
         for f in fieldnames(Agrocosm.CFTParameters))...)
end

@testset "the loss is linear in the excess millimetres" begin
    loss(mm, thr, growing, rate) =
        Agrocosm.excess_water_loss(TW(mm), TW(thr), growing, TW(rate))
    @test loss(30.0, 20.0, true, 0.01) ≈ TW(0.1)
    @test loss(40.0, 20.0, true, 0.01) ≈ TW(0.2)      # linear, not a count
    @test loss(20.0, 20.0, true, 0.01) == zero(TW)    # at the threshold, nothing
    @test loss(5.0, 20.0, true, 0.01) == zero(TW)     # below it, nothing
    @test loss(60.0, 20.0, false, 0.01) == zero(TW)   # not growing, nothing
    # A 60 mm day is not a 20.1 mm day, which is the whole reason it is not a
    # day count - and a count is not differentiable either.
    @test loss(60.0, 20.0, true, 0.01) > loss(20.1, 20.0, true, 0.01)
end

@testset "S1/S2: rate zero is inert, and a rate moves it" begin
    crop = test_model_state(init_crop(4, identity))
    state = Agrocosm.crop_prognostic(crop)
    state.phenology.harvest_recovery_fraction .= one(TW)
    state.phenology.is_growing .= Int32(1)
    Agrocosm.excess_water!(wet_cft(rate = 0.0), crop, fill(TW(80.0), 4))
    @test all(state.phenology.harvest_recovery_fraction .== one(TW))
    Agrocosm.excess_water!(wet_cft(rate = 0.01), crop, fill(TW(80.0), 4))
    @test all(state.phenology.harvest_recovery_fraction .< one(TW))
end

@testset "the damage accumulates and is irreversible" begin
    crop = test_model_state(init_crop(1, identity))
    state = Agrocosm.crop_prognostic(crop)
    state.phenology.harvest_recovery_fraction .= one(TW)
    state.phenology.is_growing .= Int32(1)
    cft = wet_cft(rate = 0.01)
    Agrocosm.excess_water!(cft, crop, fill(TW(40.0), 1))
    after_one = state.phenology.harvest_recovery_fraction[1]
    Agrocosm.excess_water!(cft, crop, fill(TW(40.0), 1))
    @test state.phenology.harvest_recovery_fraction[1] < after_one
    # A dry day afterwards does not put the crop back on its feet.
    Agrocosm.excess_water!(cft, crop, fill(TW(0.0), 1))
    @test state.phenology.harvest_recovery_fraction[1] ≈ TW(1) - TW(0.4)
    # Clamped at zero however wet it gets.
    for _ in 1:50
        Agrocosm.excess_water!(cft, crop, fill(TW(500.0), 1))
    end
    @test state.phenology.harvest_recovery_fraction[1] == zero(TW)
end

@testset "sowing restores the recovery fraction" begin
    # Without the reset, one wet season would follow the field forever - the same
    # reason `grain_set_fraction` is reset in `cultivate!`.
    source = read(joinpath(@__DIR__, "..", "..", "..", "src", "processes", "crop",
                           "cultivate.jl"), String)
    @test occursin("harvest_recovery_fraction[cell] = one(T)", source)
    @test :harvest_recovery_fraction in
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
    for cft in (Agrocosm.cft1, Agrocosm.cft2, Agrocosm.cft3, Agrocosm.cft9)
        @test cft.heavy_rain_rate == TW(0.002)
    end
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
