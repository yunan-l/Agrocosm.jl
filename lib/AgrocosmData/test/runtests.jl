using AgrocosmData
using Test

include("fixtures/fixture_data.jl")
using .FixtureData

@testset "AgrocosmData" begin
    mktempdir() do directory
        paths = FixtureData.write_fixture(directory)
        catalog = load_catalog(paths.catalog_path)

        @test dataset(catalog, :grid).path == paths.grid_path
        @test dataset(catalog, :landuse).management_bands.rainfed == Int32[1, 2]
        @test dataset(catalog, :landuse).management_bands.irrigated == Int32[3, 4]
        @test pft_index(catalog.pfts, 20) == 2
        @test pft_index(catalog.pfts, "crop_b") == 2
        @test pft_name(catalog.pfts, 10) == "crop_a"

        grid = read_grid(dataset(catalog, :grid))
        @test grid.cell_ids == Int32[0, 1, 2, 3]
        @test grid.longitude_indices == Int32[1, 3, 1, 2]
        @test grid.latitude_indices == Int32[2, 1, 1, 2]

        spatial = Float32[12 10; -1 13; 11 -1]
        compact = compact_spatial(spatial, grid, 1, 2)
        @test compact == Float32[10, 11, 12, 13]
        @test expand_to_grid(compact, grid; fill_value = -1.0f0) == spatial

        subset = select_cells(grid, [1, 3])
        @test subset.cell_ids == Int32[0, 2]
        @test compact_spatial(spatial, grid, 1, 2; selection = subset) == Float32[10, 12]

        climate_reader = climate_blocks(
            catalog, grid; co2_path = paths.co2_path, block_days = 2,
        )
        @test length(climate_reader) == 3
        climate = [block for block in climate_reader]
        @test size.(getfield.(climate, :temperature)) == fill((2, 4), 3)
        @test climate[1].co2 == Float32[369.5, 369.5]
        @test climate[2].co2 == Float32[369.5, 371.0]
        @test climate_forcing(climate[1]).co2_daily
        @test !hasproperty(climate_forcing(climate[1]), :provenance)

        eager_temp = read_compact_variable(
            dataset(catalog, :temp), grid; order = (:time, :cell), T = Float32,
        )
        @test reduce(vcat, getfield.(climate, :temperature)) == eager_temp.values
        @test reduce(vcat, getfield.(climate, :precipitation)) ==
            read_compact_variable(
                dataset(catalog, :prec), grid; order = (:time, :cell), T = Float32,
            ).values
        year_2001 = climate_blocks(
            catalog, grid; co2_path = paths.co2_path, start_year = 2001,
            end_year = 2001, block_days = 31,
        )
        @test length(year_2001) == 1
        @test read_climate_block(year_2001, 1).co2 == fill(371.0f0, 3)

        landuse = read_management(catalog, :landuse, grid, 20; years = 2000:2001)
        @test size(landuse.values) == (2, 4)
        @test landuse.values == Float32[0 0.2 0 0.3; 0.4 0 0 0.1]
        irrigated_landuse = read_management(catalog, :landuse, grid, 20; years = 2000:2001, irrigated = true)
        @test irrigated_landuse.values == landuse.values .+ 0.1f0
        crop_mask = build_crop_mask(grid, landuse.values)
        @test crop_mask.selection.cell_ids == Int32[0, 1, 3]
        @test crop_mask.fraction == Float32[0 0.2 0.3; 0.4 0 0.1]
        @test crop_mask.active == Bool[0 1 1; 1 0 1]

        residue = read_management(
            catalog,
            :residue_fraction,
            grid,
            20;
            years = 2000:2001,
            selection = crop_mask.selection,
        )
        @test size(residue.values) == (2, 3)
        @test all(0 .<= residue.values .<= 1)

        single_year = (; (
            name => read_management(
                catalog,
                name,
                grid,
                20;
                years = [2000],
                selection = crop_mask.selection,
                active = crop_mask.active[1:1, :],
            ) for name in (:sowing_date, :phu, :manure, :fertilizer, :residue_fraction)
        )...)
        crop = crop_inputs(; single_year...)
        @test crop.sdate == Int32[0, 100, 100]
        @test crop.phu == Float32[0, 1200, 1200]
        @test crop.fertilizer == Float32[0, 4, 6]
        automatic_crop = crop_inputs(
            sowing_date = single_year.sowing_date,
            phu = single_year.phu,
            residue_fraction = single_year.residue_fraction,
            fertilizer_mode = :auto,
            manure_enabled = false,
        )
        @test automatic_crop.fertilizer == zeros(Float32, 3)
        @test automatic_crop.manure == zeros(Float32, 3)

        soil = read_soil_data(catalog, grid)
        @test soil.soilcode == Int32[1, 6, 9, 14]
        @test soil.ph == Float32[6, 7, 8, 9]
        @test soil.sand == Float32[0.22, 0.58, 0.58, 0.99]
        @test size(soil.saturation) == (5, 4)
        @test soilparams(soil).soilph === soil.ph

        baseline_selection = CellSelection(1:10, 0:9)
        baseline_codes = Int32[6, 7, 9, 9, 9, 9, 9, 9, 9, 9]
        baseline_ph = Float32[6.5, 7, 7, 7, 5.5, 5.5, 5.5, 5.5, 7, 5.5]
        baseline = soil_data_from_values(baseline_codes, baseline_ph, baseline_selection)
        @test baseline.sand == Float32[0.58, 0.43, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58, 0.58]
        @test baseline.saturation[:, 1] == fill(0.404f0, 5)

        @test_throws ArgumentError validate_management(:landuse, Float32[1.1 0.0])
        @test_throws ArgumentError validate_management(:sowing_date, Float32[0 100]; active = Bool[1 1])
        @test_throws DimensionMismatch read_management(
            DatasetSpec(paths.management_path, "landfrac"; units = "%"),
            :landuse,
            grid,
            catalog.pfts,
            20,
        )
    end

    example = load_catalog(joinpath(@__DIR__, "..", "config", "catalog.example.toml"))
    @test example.pfts.ids == Int32.(1:12)
    @test pft_name(example.pfts, 4) == "tropical cereals"
    @test pft_name(example.pfts, 9) == "oil crops soybean"
    @test dataset(example, :landuse).management_bands.irrigated == Int32.(17:28)
    @test dataset(example, :sowing_date).management_bands.irrigated == Int32.(13:24)
    @test dataset(example, :residue_fraction).management_bands.irrigated == Int32.(1:12)
end
