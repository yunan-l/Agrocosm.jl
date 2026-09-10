using Agrocosm
using Test

# The cfts runner includes the wheat runner, so this brings both
# `resolved_warmup_years` and `write_batch_config` - the two halves of the
# coupling under test - into scope from one include.
#
# Idempotent: three test files need these runner scripts, and each `include`
# into Main redefines the runners' `const`s, which Julia warns may "cause
# incorrect answers". Guarding on a name each script defines makes the include
# order irrelevant.
isdefined(Main, :write_batch_config) ||
    include(joinpath(@__DIR__, "..", "..", "scripts", "run_global_cfts_cpu.jl"))

# This project runs ONE warm-up: the 600-year target-constrained calibration.
# Production continues from its pool allocation and does not spin up.
#
# The reason these are tests and not comments: every way of getting this wrong
# produces a run that completes and writes plausible output. Saying nothing in
# the config does not mean "no warm-up" - the year count falls back to 10, and
# the production phase then runs ten FREE years away from the calibration it was
# just handed, with the number appearing nowhere in the config. Observed on the
# 40-cell pre-check of 2026-09-10, whose log reported
# `warm-up completed: years=10` per rank.

@testset "Zero warm-up years is the production phase continuing from calibration" begin
    resolved = resolved_warmup_years(
        Dict("warmup_minimum_years" => 0, "warmup_maximum_years" => 0);
        target_constrained = false, has_pool_allocation = true,
    )
    @test resolved.disabled
    @test resolved.minimum_years == 0
    @test resolved.maximum_years == 0
end

@testset "Saying nothing still gives the ten-year default" begin
    # Pinned deliberately: this is the behaviour the campaign hit, and a future
    # change of the default should show up here rather than in a server log.
    resolved = resolved_warmup_years(
        Dict{String, Any}(); target_constrained = false, has_pool_allocation = true,
    )
    @test !resolved.disabled
    @test resolved.minimum_years == 10
    @test resolved.maximum_years == 10
end

@testset "The 600-year calibration is unaffected" begin
    resolved = resolved_warmup_years(
        Dict("warmup_minimum_years" => 600, "warmup_maximum_years" => 600);
        target_constrained = true, has_pool_allocation = false,
    )
    @test !resolved.disabled
    @test resolved.minimum_years == 600
    @test resolved.maximum_years == 600
end

@testset "A skip that would leave the pools uncalibrated is refused" begin
    # The calibration itself: skipping it is what makes the pools calibrated, so
    # there would be nothing to continue from.
    @test_throws ErrorException resolved_warmup_years(
        Dict("warmup_minimum_years" => 0, "warmup_maximum_years" => 0);
        target_constrained = true, has_pool_allocation = true,
    )
    # No allocation to continue from: the state would be the raw HWSD initial
    # condition, never spun up.
    @test_throws ErrorException resolved_warmup_years(
        Dict("warmup_minimum_years" => 0, "warmup_maximum_years" => 0);
        target_constrained = false, has_pool_allocation = false,
    )
    # A half-disabled reading. maximum = 0 with minimum > 0 is not a shorter
    # warm-up; it is a config whose author meant one of two different things.
    @test_throws ErrorException resolved_warmup_years(
        Dict("warmup_minimum_years" => 150, "warmup_maximum_years" => 0);
        target_constrained = false, has_pool_allocation = true,
    )
end

@testset "An absent [free_warmup] leaves production inheriting the calibration years" begin
    # `write_batch_config` overrides the warm-up keys only when a [free_warmup]
    # table is supplied, so without one the production phase carries whatever
    # [run] holds - as a FREE warm-up, because target_constrained is false for
    # production. With [run] silent that is the ten-year default; with [run]
    # carrying 600 it is 600 free years.
    base = Dict{String, Any}(
        "run" => Dict{String, Any}(
            "warmup_minimum_years" => 600, "warmup_maximum_years" => 600,
        ),
        "paths" => Dict{String, Any}("output_directory" => mktempdir()),
    )
    directory = mktempdir()
    written = write_batch_config(
        joinpath(directory, "production.toml"), base;
        output_directory = directory,
        pool_allocation = joinpath(directory, "allocation.nc"),
        production = true, target_constrained = false,
        warmup_options = nothing,
    )
    inherited = TOML.parsefile(written)["run"]
    @test inherited["warmup_maximum_years"] == 600
    @test inherited["warmup_target_constrained"] === false

    # With the table, production does not warm up at all.
    written = write_batch_config(
        joinpath(directory, "production_disabled.toml"), base;
        output_directory = directory,
        pool_allocation = joinpath(directory, "allocation.nc"),
        production = true, target_constrained = false,
        warmup_options = Dict{String, Any}(
            "warmup_minimum_years" => 0, "warmup_maximum_years" => 0,
        ),
    )
    disabled = TOML.parsefile(written)["run"]
    @test disabled["warmup_minimum_years"] == 0
    @test disabled["warmup_maximum_years"] == 0
end
