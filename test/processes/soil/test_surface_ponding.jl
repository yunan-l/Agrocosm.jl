using Test
using Agrocosm

# Ponded surface water: the store LPJmL does not have.
#
# Rejected rainfall leaves as surface runoff on the day it arrives, which is why
# this lineage cannot waterlog - measured on five real cells, root-zone
# water-filled pore space never reaches 0.90, so the aeration-stress band the
# rest of the field uses (0.85-0.95) is unreachable. `docs/32` has the table.
#
# Two properties are worth asserting and nothing else is: that capacity zero is
# bitwise the old behaviour, and that the split conserves.

const TP = Float32

@testset "capacity zero is bitwise the LPJmL behaviour" begin
    # The ablation contract. `min(rejected, 0 - 0)` is zero, so everything runs
    # off and no arithmetic differs from the version without the store.
    for rejected in (TP(0), TP(0.001), TP(12), TP(1e4))
        retained, runoff = Agrocosm.pond_rejected_water(rejected, zero(TP), zero(TP))
        @test retained === zero(TP)
        @test runoff === rejected
    end
    @test Agrocosm.ModelParameters(TP).lpjml.ponding_capacity == zero(TP)
end

@testset "the split conserves, always" begin
    # The pond is a DELAY, never a source: whatever the surface refuses is either
    # held or runs off, and the two must sum to it exactly. Asserted across the
    # corners rather than argued, because a leak here would show up as a water
    # balance residual on a global run and nowhere earlier.
    for rejected in (TP(0), TP(0.5), TP(7), TP(250)),
        ponded in (TP(0), TP(2), TP(30)),
        capacity in (TP(0), TP(1), TP(30), TP(500))
        retained, runoff = Agrocosm.pond_rejected_water(rejected, ponded, capacity)
        @test retained + runoff === rejected
        @test retained >= zero(TP)
        @test runoff >= zero(TP)
        @test ponded + retained <= max(capacity, ponded)
    end
end

@testset "a full pond passes everything on" begin
    @test Agrocosm.pond_rejected_water(TP(12), TP(5), TP(5)) === (zero(TP), TP(12))
    # And an OVERFULL pond does not drain backwards into the runoff term.
    @test Agrocosm.pond_rejected_water(TP(12), TP(9), TP(5)) === (zero(TP), TP(12))
end

@testset "the state is prognostic and reaches the water balance" begin
    # It is carried between days, so it has to be in the prognostic view or the
    # daily driver will not see yesterday's pond; and it has to be in the balance
    # or the day a pond fills reads as a leak of exactly that size.
    @test :ponding in fieldnames(Agrocosm.SoilWater)
    balance = fieldnames(Agrocosm.WaterBalance)
    @test :pond_storage_before in balance
    @test :pond_storage_after in balance
    source = read(joinpath(@__DIR__, "..", "..", "..", "src", "diagnostics",
                           "water_balance.jl"), String)
    residual = source[findfirst("water_balance.residual", source)[1]:end]
    @test occursin("pond_storage_before", residual)
    @test occursin("pond_storage_after", residual)
end
