include(joinpath(@__DIR__, "..", "..", "examples", "scripts", "run_global_cfts_cpu.jl"))

@testset "multi-CFT patch selection, state isolation, and batch configuration" begin
    high_throughput_closure = percolation_energy_closure(
        Float32[10], Float32[2.0e6], Float32[0], Float32[0], Float32[0];
        absolute_tolerance = 5.0f0, relative_tolerance = 5.0f-6,
    )
    @test high_throughput_closure.closes
    @test high_throughput_closure.maximum_relative_residual ≈ 5.0f-6
    @test !percolation_energy_closure(
        Float32[20], Float32[2.0e6], Float32[0], Float32[0], Float32[0];
        absolute_tolerance = 5.0f0, relative_tolerance = 5.0f-6,
    ).closes

    @test requested_crop_systems(Dict{String, Any}()) == [(1, false)]
    @test requested_crop_systems(Dict("cfts" => Dict(
        "cft_ids" => "all", "water_systems" => ["rainfed", "irrigated"],
    ))) == [(cft_id, irrigated) for cft_id in 1:length(CFTS) for irrigated in (false, true)]
    patch_domain = combine_patch_domains([
        PatchDomain([1], [3], [42], [1], [false], Float32[0.2]),
        PatchDomain([1], [3], [42], [2], [true], Float32[0.3]),
    ])
    @test patch_domain.cell_ids == Int32[42, 42]
    @test patch_domain.cft_ids == Int32[1, 2]
    @test patch_domain.irrigated == BitVector([false, true])

    single_catalog = catalog_from_config(Dict("paths" => Dict("input_directory" => "/tmp/input")))
    @test single_catalog.cfts.ids == Int32.(1:12)
    @test endswith(dataset(single_catalog, :landuse).path, "landuse_24cfts_2015.nc")
    @test cft_index(single_catalog.cfts, 1) == 1

    days = 30
    rainfed = initialize_simulation(
        cft1, lifecycle_initial_data(Float32);
        days, diagnostics = false, irrigation = false, fertilizer = :yes,
    )
    irrigated = initialize_simulation(
        cft1, lifecycle_initial_data(Float32);
        days, diagnostics = false, irrigation = true, fertilizer = :yes,
    )
    @test rainfed.state.prognostic.soil.water.storage !==
        irrigated.state.prognostic.soil.water.storage
    initial_irrigated_water = copy(irrigated.state.prognostic.soil.water.storage)
    run_simulation!(rainfed, lifecycle_climate(Float32, days); spinup = false)
    @test irrigated.state.prognostic.soil.water.storage == initial_irrigated_water
    run_simulation!(irrigated, lifecycle_climate(Float32, days); spinup = false)

    mktempdir() do directory
        batch_config = Dict{String, Any}(
            "paths" => Dict{String, Any}("output_directory" => directory),
            "run" => Dict{String, Any}(
                "warmup_minimum_years" => 1,
                "warmup_maximum_years" => 2,
                "require_warmup_convergence" => true,
            ),
        )
        calibration_config = write_batch_config(
            joinpath(directory, "calibration.toml"), batch_config;
            output_directory = joinpath(directory, "calibration"),
            production = false,
            target_constrained = true,
        )
        calibration_run = TOML.parsefile(calibration_config)["run"]
        @test calibration_run["warmup_minimum_years"] == 600
        @test calibration_run["warmup_maximum_years"] == 600
        @test !calibration_run["require_warmup_convergence"]

        free_config = write_batch_config(
            joinpath(directory, "production.toml"), batch_config;
            output_directory = joinpath(directory, "production"),
            pool_allocation = joinpath(directory, "allocation.nc"),
            production = true,
            target_constrained = false,
            warmup_options = Dict{String, Any}(
                "warmup_minimum_years" => 600,
                "warmup_maximum_years" => 1500,
                "require_warmup_convergence" => true,
            ),
        )
        free_run = TOML.parsefile(free_config)["run"]
        @test free_run["warmup_minimum_years"] == 600
        @test free_run["warmup_maximum_years"] == 1500
        @test free_run["require_warmup_convergence"]

        allocation_path = joinpath(directory, "allocation.nc")
        production_output = joinpath(directory, "production.nc")
        write(allocation_path, "allocation")
        write(production_output, "production")
        status_path = joinpath(directory, "batch_status.toml")
        fingerprint = _file_fingerprint(free_config)
        _write_batch_status(
            status_path;
            status = "calibration_completed",
            name = "cft_01_rainfed",
            cft_id = 1,
            irrigated = false,
            calibration_config,
            production_config = free_config,
            production_output,
            calibration_output_directory = joinpath(directory, "calibration"),
            production_output_directory = joinpath(directory, "production"),
            allocation_path,
            production_config_fingerprint = fingerprint,
            allocation_fingerprint = _file_fingerprint(allocation_path),
            calibration_completed = true,
        )
        @test _resumable_calibration(status_path, fingerprint) !== nothing
        @test _resumable_batch(status_path, fingerprint) === nothing

        _write_batch_status(
            status_path;
            status = "failed",
            name = "cft_01_rainfed",
            cft_id = 1,
            irrigated = false,
            calibration_config,
            production_config = free_config,
            production_output,
            calibration_output_directory = joinpath(directory, "calibration"),
            production_output_directory = joinpath(directory, "production"),
            allocation_path,
            production_config_fingerprint = fingerprint,
            allocation_fingerprint = _file_fingerprint(allocation_path),
            calibration_completed = true,
        )
        @test _resumable_calibration(status_path, fingerprint) !== nothing

        _write_batch_status(
            status_path;
            status = "completed",
            name = "cft_01_rainfed",
            cft_id = 1,
            irrigated = false,
            calibration_config,
            production_config = free_config,
            production_output,
            calibration_output_directory = joinpath(directory, "calibration"),
            production_output_directory = joinpath(directory, "production"),
            allocation_path,
            production_config_fingerprint = fingerprint,
        )
        @test _resumable_batch(status_path, fingerprint) !== nothing
        @test _resumable_batch(status_path, "different-config") === nothing
    end

end
