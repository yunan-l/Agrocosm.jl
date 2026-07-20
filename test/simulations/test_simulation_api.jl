using Agrocosm
using JLD2
using Test

function simulation_api_fixture(::Type{T}) where {T <: AbstractFloat}
    cells = 1
    layers = 5
    initial = (
        coords = [1],
        latitude = T[45],
        crop = (
            sdate = Int32[1],
            phu = T[543],
            manure = zeros(T, cells),
            fertilizer = T[24.55],
            residuefrac = T[0.67],
        ),
        soilparam = (
            soilph = T[6.5],
            w_sat = fill(T(0.45), layers, cells),
            sand = T[0.4],
            clay = T[0.2],
            tdiff_0 = T[0.7],
            tdiff_15 = T[0.75],
            soildepth = T[200, 300, 500, 1000, 1000],
        ),
        initialLPJmL = (u0 = (
            swc = reshape(T[57.41, 55.32, 126.13, 274.59, 285.71], layers, cells),
            litc = reshape(T[0.13, 187.5, 225.36], 3, cells),
            fastc = fill(T(10), layers, cells),
            slowc = fill(T(100), layers, cells),
            litn = reshape(T[0.0047, 6.47, 9.47], 3, cells),
            fastn = fill(T(1), layers, cells),
            slown = fill(T(10), layers, cells),
        ),),
    )
    climate = (
        temp_spinup = fill(T(10), 365, cells),
        temp = fill(T(15), 3, cells),
        prec = fill(T(1), 3, cells),
        swdown = fill(T(180), 3, cells),
        lwnet = fill(T(-40), 3, cells),
        windspeed = fill(T(2), 3, cells),
        co2 = T[400],
    )
    return initial, climate
end

function climate_block(::Type{T}, days, temperature, precipitation) where {T <: AbstractFloat}
    return (
        temp_spinup = fill(T(10), 365, 1),
        temp = fill(T(temperature), days, 1),
        prec = fill(T(precipitation), days, 1),
        swdown = fill(T(180), days, 1),
        lwnet = fill(T(-40), days, 1),
        windspeed = fill(T(2), days, 1),
        co2 = T[400],
    )
end

@testset "Annual CO₂ forcing length is validated before kernel launch" begin
    initial, _ = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float32, days = 366, auto_fertilizer = false,
    )
    incomplete = climate_block(Float32, 366, 15, 1)
    @test_throws DimensionMismatch run_simulation!(
        simulation, incomplete; spinup = false,
    )
end

@testset "High-level crop simulation API" begin
    initial, climate = simulation_api_fixture(Float32)
    simulation = initialize_simulation(
        cft1, initial;
        indices = [1],
        T = Float64,
        device = identity,
        days = 3,
        auto_fertilizer = false,
    )

    @test simulation.crop === simulation.state.crop
    @test simulation.output === simulation.state.output
    @test eltype(simulation.crop.canopy.lai) == Float64
    @test eltype(simulation.water_balance.residual) == Float64

    returned = run_simulation!(simulation, climate; spinup = false)
    @test returned === simulation
    @test simulation.simulated_days == 3
    @test size(simulation.output.crop.npp) == (3, 1)
    @test all(isfinite, simulation.output.crop.npp)

    summary = simulation_summary(simulation)
    @test summary.precision == Float64
    @test summary.cells == 1
    @test summary.simulated_days == 3
    @test isfinite(summary.crop.cumulative_npp)
    @test_throws ArgumentError run_simulation!(simulation, climate; spinup = false)
end

@testset "Multiple climate blocks preserve the continuous daily timeline" begin
    initial, _ = simulation_api_fixture(Float32)
    first_block = climate_block(Float32, 2, 15, 1)
    second_block = climate_block(Float32, 2, 17, 3)
    continuous = (
        temp_spinup = first_block.temp_spinup,
        temp = vcat(first_block.temp, second_block.temp),
        prec = vcat(first_block.prec, second_block.prec),
        swdown = vcat(first_block.swdown, second_block.swdown),
        lwnet = vcat(first_block.lwnet, second_block.lwnet),
        windspeed = vcat(first_block.windspeed, second_block.windspeed),
        co2 = Float32[400],
    )

    create() = initialize_simulation(
        cft1, initial;
        indices = [1], T = Float64, days = 4, auto_fertilizer = false,
    )
    chunked = create()
    single = create()
    run_simulation!(chunked, [first_block, second_block]; spinup = false)
    run_simulation!(single, continuous; spinup = false)

    @test chunked.simulated_days == 4
    @test chunked.output.crop.npp ≈ single.output.crop.npp
    @test chunked.output.calendar.sowing_callback == single.output.calendar.sowing_callback
    @test findall(!iszero, vec(chunked.output.calendar.sowing_callback)) == [1]
    @test chunked.crop.carbon.organs ≈ single.crop.carbon.organs
    @test chunked.soil.water.storage ≈ single.soil.water.storage
    @test chunked.water_balance.precipitation == reshape(Float64[1, 1, 3, 3], 4, 1)

    mktempdir() do directory
        first_path = joinpath(directory, "climate_1.jld2")
        second_path = joinpath(directory, "climate_2.jld2")
        JLD2.jldsave(first_path; climate = first_block)
        JLD2.jldsave(second_path; climate = second_block)
        from_files = create()
        run_simulation!(from_files, [first_path, second_path]; spinup = false)
        @test from_files.output.crop.npp ≈ single.output.crop.npp
        @test from_files.soil.water.storage ≈ single.soil.water.storage
    end
end
