using Agrocosm
using Terrarium
using Test
import RingGrids

# `surface_climate_inputs` packs tabulated climate series into Terrarium time-varying input sources
# (Oceananigans `FieldTimeSeries` + `FieldTimeSeriesInputSource`). These check that the forcing reaches
# the model state, that per-column (matrix) series are accepted on a global grid, and that shape
# mismatches error.
@testset "Surface climate input sources" begin
    ndays = 4
    temp = [5.0, 10.0, 15.0, 20.0]                                   # °C, one per day
    times = (0:ndays - 1) .* Terrarium.seconds_per_day(Float64)      # 0, 1 day, 2 days, 3 days

    @testset "vector series drives the model input" begin
        grid = ColumnGrid(CPU(), UniformSpacing(Δz = 0.1, N = 5))
        model = CropModel(grid, crop_pft("maize"))
        inputs = surface_climate_inputs(grid, times; air_temperature = temp)
        @test inputs isa Terrarium.InputSources
        integrator = initialize(model; inputs, initializers = (temperature = 10.0,))
        # at t = 0 the interpolated input equals the first sample
        @test interior(integrator.state.air_temperature)[1, 1, 1] ≈ temp[1]
    end

    @testset "length mismatch errors" begin
        grid = ColumnGrid(CPU(), UniformSpacing(Δz = 0.1, N = 5))
        @test_throws ArgumentError surface_climate_inputs(grid, times; air_temperature = [1.0, 2.0])
    end

    @testset "per-column (matrix) series on a global grid" begin
        rings = RingGrids.FullGaussianGrid(1)
        grid = ColumnRingGrid(CPU(), Float64, UniformSpacing(Δz = 0.1, N = 5), rings)
        ncolumns = length(RingGrids.get_lonlats(grid.rings)[1])
        matrix = Float64[10 + c for t in 1:ndays, c in 1:ncolumns]   # column c constant at 10 + c
        inputs = surface_climate_inputs(grid, times; air_temperature = matrix)
        @test inputs isa Terrarium.InputSources
        # wrong number of columns is rejected
        @test_throws ArgumentError surface_climate_inputs(grid, times; air_temperature = matrix[:, 1:1])
    end
end
