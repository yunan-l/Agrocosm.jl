using NCDatasets

include(joinpath(@__DIR__, "..", "..", "examples", "scripts", "run_global_cfts_cpu.jl"))

@testset "multi-patch state isolation and yield bands" begin
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
        batch_paths = String[]
        for (index, values) in enumerate((Float32[1 2], Float32[3 4]))
            path = joinpath(directory, "batch_$index.nc")
            NCDataset(path, "c") do dataset
                defDim(dataset, "longitude", 2)
                defDim(dataset, "latitude", 1)
                defDim(dataset, "time", 1)
                defVar(dataset, "longitude", Float32, ("longitude",))[:] = Float32[0, 1]
                defVar(dataset, "latitude", Float32, ("latitude",))[:] = Float32[0]
                defVar(dataset, "time", Int32, ("time",))[:] = Int32[2015]
                defVar(dataset, "crop_yield", Float32, ("longitude", "latitude", "time"))[:, :, :] =
                    reshape(values, 2, 1, 1)
            end
            push!(batch_paths, path)
        end
        output = write_cft_yield(
            joinpath(directory, "yield.nc"), batch_paths, [(1, false), (1, true)],
        )
        NCDataset(output, "r") do dataset
            @test size(dataset["yield"]) == (2, 1, 2, 1)
            @test dataset["yield"][:, :, 1, :] == reshape(Float32[1, 2], 2, 1, 1)
            @test dataset["yield"][:, :, 2, :] == reshape(Float32[3, 4], 2, 1, 1)
            @test dataset["pft_id"][:] == Int32[1, 1]
            @test dataset["irrigated"][:] == Int8[0, 1]
        end
    end
end
