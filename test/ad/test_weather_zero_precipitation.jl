using Agrocosm
using Enzyme
using Test

include(joinpath(@__DIR__, "..", "helpers", "weather_attribution_fixture.jl"))

for cft_id in (isempty(ARGS) ? (1, 3) : Tuple(parse.(Int, ARGS)))
    @testset "CFT $cft_id weather AD with dry days" begin
        case = weather_attribution_fixture(cft_id; T = Float32)
        (; forcing, state, cft, parameters, climate, days, harvest_day) = case
        dry_days = collect(last(days):-2:first(days))
        forcing[dry_days, 1, 2] .= 0.0f0
        original = copy(forcing)
        reference = weather_harvest_replay(forcing, state, cft, parameters,
            climate, days, harvest_day)
        @test reference.schedule_matches && !reference.failed
        @test reference.yield > 0
        @info "Dry-day weather regression" cft_id production_yield = reference.yield
        flush(stderr)

        # Both rain and melt are zero on dry days. The inactive enthalpy
        # fractions must not introduce 0/0 into the reverse weather gradient.
        result = enzyme_weather_harvest_gradient(forcing, state, cft, parameters,
            climate, days, harvest_day; block_days = 4)
        @test all(isfinite, result.gradient)
        @test all(isfinite, result.gradient[dry_days, 1, 2])
        @test result.primal ≈ reference.yield rtol = 1e-3 atol = 1e-5
        @test result.reverse_primal ≈ result.primal rtol = 1e-5 atol = 1e-6
        @test forcing == original
        @test all(iszero, result.gradient[1:first(days)-1, :, :])
        @test all(iszero, result.gradient[harvest_day:end, :, :])

        # A central difference at zero rain would cross the forcing domain.
        # Probe wet days instead, while reverse propagation still crosses dry days.
        direction = zeros(Float32, size(forcing))
        wet_days = [day for day in days if forcing[day, 1, 2] > 0]
        direction[wet_days, 1, 2] .= 1.0f0
        projection = sum(result.gradient .* direction)
        delta = 0.05f0
        plus = weather_harvest_replay(forcing .+ delta .* direction, state,
            cft, parameters, climate, days, harvest_day)
        minus = weather_harvest_replay(forcing .- delta .* direction, state,
            cft, parameters, climate, days, harvest_day)
        @test plus.schedule_matches && minus.schedule_matches
        @test projection ≈ (plus.yield - minus.yield) / (2delta) rtol = 0.03 atol = 2e-5
    end
end
