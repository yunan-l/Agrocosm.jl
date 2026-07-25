using JLD2
using NCDatasets

function stream_fixture_climate(climate, indices)
    return ClimateDataLoader(climate, indices, identity; T = Float32)
end

function split_stream_climate(climate, block_days, days = size(climate.temp, 1))
    return [(
        temp = copy(climate.temp[rows, :]),
        prec = copy(climate.prec[rows, :]),
        sw = copy(climate.sw[rows, :]),
        lw = copy(climate.lw[rows, :]),
        co2 = ndims(climate.co2) == 1 ? fill(climate.co2[1], length(rows)) :
            copy(climate.co2[rows, :]),
        co2_daily = true,
    ) for first_day in 1:block_days:days
      for rows in (first_day:min(days, first_day + block_days - 1),)]
end

@testset "streamed ten-cell daily output equals in-memory output" begin
    example_dir = joinpath(@__DIR__, "..", "..", "examples")
    initial = JLD2.load(joinpath(example_dir, "initial_wheat.jld2"), "initial_data")
    raw_climate = JLD2.load(joinpath(example_dir, "climate_2000_2009.jld2"), "climate")
    cells = length(initial.latitude)
    indices = collect(1:cells)
    days = 365
    climate = stream_fixture_climate(raw_climate, indices)

    baseline = initialize_simulation(
        cft1, initial;
        indices, T = Float32, days, diagnostics = false, fertilizer = :yes,
    )
    run_simulation!(baseline, climate; end_day = days, spinup = false)

    chunks = OutputChunk[]
    stream = OutputStream(
        [
            OutputVariable(:crop, :gpp),
            OutputVariable(:crop, :npp),
            OutputVariable(:crop, :lai),
            OutputVariable(:crop, :yield),
        ];
        frequency = :daily,
        writer = chunk -> push!(chunks, chunk),
        cell_ids = indices,
    )
    streamed = initialize_simulation(
        cft1, initial;
        indices, T = Float32, days, diagnostics = false, fertilizer = :yes,
    )
    blocks = split_stream_climate(climate, 73, days)
    run_simulation!(streamed, blocks; spinup = false, output_stream = stream)

    daily_chunks = filter(chunk -> chunk.frequency === :daily, chunks)
    annual_chunks = filter(chunk -> chunk.frequency === :annual, chunks)
    @test reduce(vcat, [chunk.values[:crop_gpp] for chunk in daily_chunks]) ==
        Array(baseline.output.crop.gpp)
    @test reduce(vcat, [chunk.values[:crop_npp] for chunk in daily_chunks]) ==
        Array(baseline.output.crop.npp)
    @test reduce(vcat, [chunk.values[:crop_lai] for chunk in daily_chunks]) ==
        Array(baseline.output.crop.lai)
    @test reduce(vcat, [chunk.values[:crop_yield] for chunk in annual_chunks]) ==
        Array(baseline.output.crop.yield)
    @test all(length(chunk.time) <= 73 for chunk in daily_chunks)
    @test all(Set(keys(chunk.values)) == Set((:crop_gpp, :crop_npp, :crop_lai))
              for chunk in daily_chunks)
    @test Agrocosm._output_timeseries_empty(streamed.output)

    diagnostic_simulation = initialize_simulation(
        cft1, initial;
        indices, T = Float32, days = 1, diagnostics = true, fertilizer = :no,
    )
    diagnostic_stream = OutputStream(
        [OutputVariable(:crop, :npp)]; cell_ids = indices,
    )
    @test_throws ArgumentError run_simulation!(
        diagnostic_simulation, [blocks[1]]; spinup = false,
        output_stream = diagnostic_stream,
    )
end

@testset "monthly and annual stream aggregation" begin
    monthly_chunks = OutputChunk[]
    monthly = OutputStream(
        [
            OutputVariable(:crop, :npp; reduction = :sum),
            OutputVariable(:crop, :lai; reduction = :mean),
        ];
        frequency = :monthly,
        writer = chunk -> push!(monthly_chunks, chunk),
        cell_ids = 1:2,
    )
    output = init_output(Float32, 2, identity)
    Agrocosm.prepare_output_block!(output, 30, 0)
    output.crop.npp .= 2
    output.crop.lai .= 4
    consume_output!(monthly, output, 1)
    clear_output_timeseries!(output)
    Agrocosm.prepare_output_block!(output, 30, 0)
    output.crop.npp .= 2
    output.crop.lai .= 4
    consume_output!(monthly, output, 31)
    finish_output_stream!(monthly, 60)
    @test getfield.(monthly_chunks, :time) == [[31], [59], [60]]
    @test monthly_chunks[1].values[:crop_npp] == fill(62.0f0, 1, 2)
    @test monthly_chunks[1].values[:crop_lai] == fill(4.0f0, 1, 2)
    @test monthly_chunks[2].values[:crop_npp] == fill(56.0f0, 1, 2)
    @test monthly_chunks[3].values[:crop_npp] == fill(2.0f0, 1, 2)

    annual_chunks = OutputChunk[]
    annual = OutputStream(
        [
            OutputVariable(:crop, :npp; reduction = :sum),
            OutputVariable(:crop, :yield),
        ];
        frequency = :annual,
        writer = chunk -> push!(annual_chunks, chunk),
        cell_ids = 1:2,
    )
    clear_output_timeseries!(output)
    Agrocosm.prepare_output_block!(output, 365, 1)
    output.crop.npp .= 1
    output.crop.yield .= reshape(Float32[10, 20], 1, :)
    consume_output!(annual, output, 1)
    finish_output_stream!(annual, 365)
    @test length(annual_chunks) == 1
    @test annual_chunks[1].time == [365]
    @test annual_chunks[1].values[:crop_npp] == fill(365.0f0, 1, 2)
    @test annual_chunks[1].values[:crop_yield] == reshape(Float32[10, 20], 1, :)
end

@testset "stream block writers" begin
    chunk = OutputChunk(
        1, :daily, [1, 2], Int32[7, 9],
        Dict{Symbol, Any}(:crop_npp => Float32[1 2; 3 4]),
    )
    mktempdir() do directory
        jld_path = JLD2BlockWriter(directory; prefix = "jld")(chunk)
        @test isfile(jld_path)
        @test JLD2.load(jld_path, "values")[:crop_npp] == chunk.values[:crop_npp]

        nc_path = NetCDFBlockWriter(directory; prefix = "nc")(chunk)
        @test isfile(nc_path)
        NCDataset(nc_path, "r") do dataset
            @test dataset["time"][:] == Int32[1, 2]
            @test dataset["cell_id"][:] == Int32[7, 9]
            @test dataset["crop_npp"][:, :] == chunk.values[:crop_npp]
        end
    end
end
