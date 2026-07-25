@testset "HWSD C/N preprocessing" begin
    organic_carbon = fill(1.0f0, 7, 2)
    total_nitrogen = fill(1.0f0, 7, 2)
    bulk_density = fill(1.0f0, 7, 2)
    coarse_fragments = fill(0.0f0, 7, 2)
    stocks = hwsd_layer_stocks(
        organic_carbon, total_nitrogen, bulk_density, coarse_fragments,
    )
    @test stocks.carbon[1, 1] == 2000.0f0
    @test stocks.carbon[6, 1] == 5000.0f0
    @test stocks.nitrogen[1, 1] == 200.0f0
    @test stocks.nitrogen[6, 1] == 500.0f0

    carbon_layers = remap_hwsd_layers(stocks.carbon)
    @test carbon_layers.values[:, 1] == Float32[2000, 3000, 5000, 10000, 10000]
    @test sum(carbon_layers.values[1:4, 1]) == sum(stocks.carbon[:, 1])
    @test carbon_layers.uncertain[:, 1] == Bool[0, 0, 0, 0, 1]
    missing_deep = remap_hwsd_layers(stocks.carbon; deep_rule = :missing)
    @test isnan(missing_deep.values[5, 1])
    @test missing_deep.uncertain[5, 1]

    organic_carbon[:, 2] .= 2
    selection = CellSelection(1:1, Int32[42])
    targets = preprocess_hwsd_cn(
        organic_carbon,
        total_nitrogen,
        bulk_density,
        coarse_fragments,
        Int32[1, 1],
        Float64[1, 3],
        selection,
    )
    @test targets.soil_organic_carbon[:, 1] ==
        Float32[3500, 5250, 8750, 17500, 17500]
    @test targets.total_nitrogen[:, 1] ==
        Float32[200, 300, 500, 1000, 1000]
    @test targets.coverage == ones(Float32, 5, 1)
    @test targets.uncertain[:, 1] == Bool[0, 0, 0, 0, 1]
    @test targets.provenance.source_version == "HWSD v2.01"
    @test targets.provenance.deep_rule == :extend_deepest_density
    @test targets.provenance.conservation.carbon_g ==
        vec(sum(targets.soil_organic_carbon .* 4; dims = 2))

    grid = GridIndex(
        Float64[-0.25, 0.25], Float64[0.25, 0.75],
        Int32[0 2; 1 3], Int32[0, 1, 2, 3],
        Int32[1, 2, 1, 2], Int32[1, 1, 2, 2],
    )
    mapping = hwsd_tile_mapping(
        Float64[-0.375, -0.125, 0.125, 0.375],
        Float64[0.125, 0.375, 0.625, 0.875],
        grid,
    )
    @test reshape(mapping.target_indices, 4, 4) == Int32[
        1 1 3 3
        1 1 3 3
        2 2 4 4
        2 2 4 4
    ]
    @test all(mapping.pixel_area .> 0)
    @test isapprox(
        sum(reshape(mapping.pixel_area, 4, 4)[1:2, 1:2]),
        AgrocosmData._EARTH_RADIUS_M^2 * deg2rad(0.5) *
            (sin(deg2rad(0.5)) - sin(deg2rad(0.0)));
        rtol = 1e-14,
    )

    missing_carbon = copy(organic_carbon)
    missing_carbon[:, 1] .= NaN
    partial = preprocess_hwsd_cn(
        missing_carbon,
        total_nitrogen,
        bulk_density,
        coarse_fragments,
        Int32[1, 1],
        Float64[1, 3],
        selection;
        minimum_coverage = 0.7,
    )
    @test partial.coverage[:, 1] == fill(0.75f0, 5)
    @test partial.soil_organic_carbon[:, 1] ==
        Float32[4000, 6000, 10000, 20000, 20000]
    @test all(partial.uncertain)

    mktempdir() do directory
        path = write_soil_cn_targets(joinpath(directory, "hwsd_cn.nc"), targets)
        restored = read_soil_cn_targets(path)
        @test restored.selection.cell_ids == targets.selection.cell_ids
        @test restored.layer_bounds == targets.layer_bounds
        @test restored.soil_organic_carbon == targets.soil_organic_carbon
        @test restored.total_nitrogen == targets.total_nitrogen
        @test restored.coverage == targets.coverage
        @test restored.uncertain == targets.uncertain
        @test restored.provenance.source_version == targets.provenance.source_version
    end

    @test_throws ArgumentError hwsd_layer_stocks(
        fill(101.0, 7, 1), fill(1.0, 7, 1), fill(1.0, 7, 1), fill(0.0, 7, 1),
    )
    @test_throws ArgumentError remap_hwsd_layers(stocks.carbon; deep_rule = :invent)
end
