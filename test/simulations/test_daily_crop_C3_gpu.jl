using Agrocosm
using CUDA
using Test

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

const C3_E2E_DAYS = 730
const C3_E2E_CELLS = 2

to_backend(::typeof(identity), values) = copy(values)
to_backend(device, values) = device(values)

function c3_e2e_initial_data(device)
    cells = C3_E2E_CELLS
    layers = 5

    crop = (
        sdate = to_backend(device, Int32[100, 105]),
        phu = to_backend(device, Float32[543, 620]),
        manure = to_backend(device, zeros(Float32, cells)),
        fertilizer = to_backend(device, Float32[24.55, 30]),
        residuefrac = to_backend(device, Float32[0.67, 0.55]),
    )
    soil_parameters = (
        ph = to_backend(device, Float32[6.5, 6.2]),
        w_sat = to_backend(device, repeat(Float32[0.404, 0.43]', layers, 1)),
        sand = to_backend(device, reshape(Float32[0.58, 0.46], 1, cells)),
        clay = to_backend(device, reshape(Float32[0.27, 0.22], 1, cells)),
        tdiff_0 = to_backend(device, Float32[0.78, 0.70]),
        tdiff_15 = to_backend(device, Float32[0.808, 0.75]),
        # Keep this host-side: init_soil transfers the immutable layer geometry.
        soildepth = Float32[200, 300, 500, 1000, 1000],
    )

    swc = Float32[57.41 62.0; 55.32 68.0; 126.13 140.0; 274.59 260.0; 285.71 275.0]
    litc = Float32[0.13 0.2; 187.5 165.0; 225.36 205.0]
    fastc = Float32[548.97 510.0; 368.27 350.0; 313.79 300.0; 377.55 360.0; 344.65 330.0]
    slowc = Float32[1218.62 1150.0; 753.33 720.0; 660.10 640.0; 792.63 760.0; 736.38 710.0]
    litn = Float32[0.0047 0.007; 6.47 5.8; 9.47 8.7]
    fastn = Float32[36.60 34.0; 24.55 23.0; 20.92 20.0; 25.17 24.0; 22.98 22.0]
    slown = Float32[81.24 77.0; 50.22 48.0; 44.01 42.5; 52.84 50.5; 49.09 47.0]
    u0 = (
        swc = to_backend(device, swc),
        litc = to_backend(device, litc),
        fastc = to_backend(device, fastc),
        slowc = to_backend(device, slowc),
        litn = to_backend(device, litn),
        fastn = to_backend(device, fastn),
        slown = to_backend(device, slown),
    )

    return (
        latitude = to_backend(device, Float32[45, 47]),
        soilparams = soil_parameters,
        ModelState = (crop = crop, u0 = u0),
    )
end

function c3_e2e_climate_host()
    days = C3_E2E_DAYS
    cells = C3_E2E_CELLS
    temperature = zeros(Float32, days, cells)
    precipitation = zeros(Float32, days, cells)
    shortwave = zeros(Float32, days, cells)
    longwave = zeros(Float32, days, cells)
    wind = zeros(Float32, days, cells)

    for cell in 1:cells, day in 1:days
        day_of_year = mod1(day, 365)
        phase = 2.0f0 * Float32(pi) * Float32(day_of_year - 109 - 5 * (cell - 1)) / 365.0f0
        seasonal = sin(phase)
        temperature[day, cell] = 9.0f0 + 11.0f0 * seasonal + 0.4f0 * (cell - 1)
        shortwave[day, cell] = max(20.0f0, 175.0f0 + 125.0f0 * seasonal)
        longwave[day, cell] = -45.0f0 + 12.0f0 * seasonal
        wind[day, cell] = 1.8f0 + 0.2f0 * (cell - 1)
        precipitation[day, cell] = if mod(day + 2 * cell, 11) == 0
            11.0f0
        elseif mod(day + cell, 5) == 0
            2.5f0
        else
            0.0f0
        end
    end

    return (
        temp = temperature,
        prec = precipitation,
        sw = shortwave,
        lw = longwave,
        wind = wind,
        co2 = Float32[370, 373],
    )
end

function c3_e2e_climate(device, host)
    return (
        temp = to_backend(device, host.temp),
        prec = to_backend(device, host.prec),
        sw = to_backend(device, host.sw),
        lw = to_backend(device, host.lw),
        wind = to_backend(device, host.wind),
        co2 = to_backend(device, host.co2),
    )
end

function initialize_c3_e2e_climbuf!(climbuf, host_climate, device)
    # Use the same deterministic one-year history on both backends without
    # adding another 365 kernel launches to this end-to-end regression.
    climbuf.atemp .= to_backend(device, host_climate.temp[1:365, :])
    climbuf.temp .= to_backend(device, host_climate.temp[335:365, :])
    climbuf.atemp_mean .= to_backend(
        device,
        vec(sum(host_climate.temp[1:365, :]; dims = 1) ./ 365.0f0),
    )
    climbuf.V_req .= 0.0f0
    return nothing
end

function run_c3_e2e(device, host_climate)
    initial_data = c3_e2e_initial_data(device)
    climbuf, crop, pet, soil, managed_land, weather, output = init_states!(
        cft1, initial_data, C3_E2E_CELLS, device,
    )
    climate = c3_e2e_climate(device, host_climate)
    initialize_c3_e2e_climbuf!(climbuf, host_climate, device)

    water = init_water_balance(C3_E2E_DAYS, C3_E2E_CELLS, device)
    nitrogen = init_nitrogen_balance(C3_E2E_DAYS, C3_E2E_CELLS, device)
    carbon = init_carbon_balance(C3_E2E_DAYS, C3_E2E_CELLS, device)
    thermal = init_thermal_balance(C3_E2E_DAYS, C3_E2E_CELLS, device)

    daily_crop_C3!(
        1, C3_E2E_DAYS,
        cft1, climate, climbuf, crop, pet, soil, managed_land, weather, output;
        irrigation = false,
        manure = false,
        auto_fertilizer = false,
        nitrogen_limit_vmax = false,
        water_balance = water,
        nitrogen_balance = nitrogen,
        carbon_balance = carbon,
        thermal_balance = thermal,
    )

    return (; climbuf, crop, pet, soil, managed_land, weather, output,
            water, nitrogen, carbon, thermal)
end

host_array(values) = Array(values)

function test_float_equivalence(gpu_values, cpu_values;
                                rtol = 1.0f-3, atol = 5.0f-4,
                                label = "unnamed")
    gpu_host = host_array(gpu_values)
    cpu_host = host_array(cpu_values)
    @test size(gpu_host) == size(cpu_host)
    @test all(isfinite, gpu_host)
    @test all(isfinite, cpu_host)
    # Array `isapprox` uses a norm, so thousands of individually acceptable
    # daily Float32 differences can fail after accumulation. Model equivalence
    # is a pointwise requirement: every day/cell must satisfy the tolerance.
    matches = isapprox.(gpu_host, cpu_host; rtol = rtol, atol = atol)
    if !all(matches)
        first_bad = findfirst(.!matches)
        gpu_value = gpu_host[first_bad]
        cpu_value = cpu_host[first_bad]
        absolute_error = abs(gpu_value - cpu_value)
        relative_error = absolute_error / max(abs(cpu_value), atol)
        @info "CPU/GPU first pointwise mismatch" label first_bad gpu_value cpu_value absolute_error relative_error rtol atol
    end
    @test all(matches)
    return nothing
end

function test_balance_equivalence(gpu_balance, cpu_balance;
                                  group, rtol, atol,
                                  field_atol = NamedTuple())
    for field in fieldnames(typeof(cpu_balance))
        comparison_atol = hasproperty(field_atol, field) ?
            getproperty(field_atol, field) : atol
        test_float_equivalence(
            getproperty(gpu_balance, field),
            getproperty(cpu_balance, field);
            rtol = rtol,
            atol = comparison_atol,
            label = "$group.$field",
        )
    end
    return nothing
end

function callback_days(values, cell)
    host = host_array(values)
    return findall(!iszero, view(host, :, cell))
end

function log_first_crop_divergence(cpu, gpu;
                                   rtol = 1.0f-5, atol = 1.0f-6)
    for cell in 1:C3_E2E_CELLS
        cpu_sowing = callback_days(cpu.output.calendar.sowing_callback, cell)
        gpu_sowing = callback_days(gpu.output.calendar.sowing_callback, cell)
        cpu_harvest = callback_days(cpu.output.calendar.harvest_callback, cell)
        gpu_harvest = callback_days(gpu.output.calendar.harvest_callback, cell)
        @info "CPU/GPU crop event days" cell cpu_sowing gpu_sowing cpu_harvest gpu_harvest
    end

    # Follow the daily causal chain from phenology/photosynthesis to crop state.
    for field in (:fphu, :lambda, :gpp, :respiration, :npp, :biomass, :lai)
        cpu_values = host_array(getproperty(cpu.output.crop, field))
        gpu_values = host_array(getproperty(gpu.output.crop, field))
        matches = isapprox.(gpu_values, cpu_values; rtol = rtol, atol = atol)
        first_bad = findfirst(.!matches)
        first_bad === nothing && continue
        gpu_value = gpu_values[first_bad]
        cpu_value = cpu_values[first_bad]
        absolute_error = abs(gpu_value - cpu_value)
        relative_error = absolute_error / max(abs(cpu_value), atol)
        @info "CPU/GPU key crop divergence" field first_bad gpu_value cpu_value absolute_error relative_error rtol atol
    end
    return nothing
end

function test_exact_equivalence(gpu_values, cpu_values; label)
    gpu_host = host_array(gpu_values)
    cpu_host = host_array(cpu_values)
    matches = gpu_host .== cpu_host
    if !all(matches)
        first_bad = findfirst(.!matches)
        gpu_value = gpu_host[first_bad]
        cpu_value = cpu_host[first_bad]
        @info "CPU/GPU first exact mismatch" label first_bad gpu_value cpu_value
    end
    @test all(matches)
    return nothing
end

@testset "CUDA C3 rainfed wheat 365/730-day end-to-end equivalence" begin
    host_climate = c3_e2e_climate_host()
    cpu = run_c3_e2e(identity, host_climate)
    gpu = run_c3_e2e(CuArray, host_climate)
    synchronize()

    @test size(cpu.output.crop.npp) == (730, C3_E2E_CELLS)
    @test size(gpu.output.crop.npp) == (730, C3_E2E_CELLS)
    @test maximum(cpu.output.crop.npp) > 0.0f0
    @test sum(cpu.output.calendar.sowing_callback) > 0
    @test sum(cpu.output.calendar.harvest_callback) > 0

    log_first_crop_divergence(cpu, gpu)

    daily_float_fields = (
        :gpp, :npp, :lambda, :potential_vmax, :vmax,
        :nitrogen_limitation, :respiration, :biomass, :lai,
        :storage_carbon, :fphu, :water_deficit,
    )
    for field in daily_float_fields
        cpu_values = getproperty(cpu.output.crop, field)
        gpu_values = getproperty(gpu.output.crop, field)
        # Explicit one-year checkpoint and complete two-year trajectory.
        test_float_equivalence(
            gpu_values[1:365, :], cpu_values[1:365, :];
            label = "output.crop.$field[1:365]",
        )
        test_float_equivalence(
            gpu_values, cpu_values; label = "output.crop.$field[1:730]",
        )
    end
    test_float_equivalence(
        gpu.output.crop.yield, cpu.output.crop.yield;
        label = "output.crop.yield",
    )

    for field in (:growing_mask,)
        test_exact_equivalence(
            getproperty(gpu.output.crop, field),
            getproperty(cpu.output.crop, field);
            label = "output.crop.$field",
        )
    end
    for field in (
        :harvesting_mask, :harvesting_year, :harvest_date,
        :sowing_callback, :harvest_callback,
    )
        test_exact_equivalence(
            getproperty(gpu.output.calendar, field),
            getproperty(cpu.output.calendar, field);
            label = "output.calendar.$field",
        )
    end

    crop_state_fields = (
        ("crop.carbon.organs", cpu.crop.carbon.organs, gpu.crop.carbon.organs),
        ("crop.nitrogen.total", cpu.crop.nitrogen.total, gpu.crop.nitrogen.total),
        ("crop.nitrogen.leaf", cpu.crop.nitrogen.leaf, gpu.crop.nitrogen.leaf),
        ("crop.nitrogen.root", cpu.crop.nitrogen.root, gpu.crop.nitrogen.root),
        ("crop.nitrogen.pool", cpu.crop.nitrogen.pool, gpu.crop.nitrogen.pool),
        ("crop.nitrogen.storage", cpu.crop.nitrogen.storage, gpu.crop.nitrogen.storage),
        ("crop.canopy.lai", cpu.crop.canopy.lai, gpu.crop.canopy.lai),
        ("crop.phenology.fphu", cpu.crop.phenology.fphu, gpu.crop.phenology.fphu),
    )
    for (label, cpu_values, gpu_values) in crop_state_fields
        test_float_equivalence(gpu_values, cpu_values; label = label)
    end

    soil_state_fields = (
        ("soil.water.storage", cpu.soil.water.storage, gpu.soil.water.storage),
        ("soil.water.ice_storage", cpu.soil.water.ice_storage, gpu.soil.water.ice_storage),
        ("soil.thermal.temperature", cpu.soil.thermal.temperature, gpu.soil.thermal.temperature),
        ("soil.thermal.enthalpy", cpu.soil.thermal.enthalpy, gpu.soil.thermal.enthalpy),
        ("soil.carbon.litter", cpu.soil.carbon.litter, gpu.soil.carbon.litter),
        ("soil.carbon.fast", cpu.soil.carbon.fast, gpu.soil.carbon.fast),
        ("soil.carbon.slow", cpu.soil.carbon.slow, gpu.soil.carbon.slow),
        ("soil.nitrogen.litter", cpu.soil.nitrogen.litter, gpu.soil.nitrogen.litter),
        ("soil.nitrogen.fast", cpu.soil.nitrogen.fast, gpu.soil.nitrogen.fast),
        ("soil.nitrogen.slow", cpu.soil.nitrogen.slow, gpu.soil.nitrogen.slow),
        ("soil.nitrogen.nitrate", cpu.soil.nitrogen.nitrate, gpu.soil.nitrogen.nitrate),
        ("soil.nitrogen.ammonium", cpu.soil.nitrogen.ammonium, gpu.soil.nitrogen.ammonium),
    )
    for (label, cpu_values, gpu_values) in soil_state_fields
        test_float_equivalence(gpu_values, cpu_values; label = label)
    end
    @test host_array(gpu.soil.thermal.initialized) ==
        host_array(cpu.soil.thermal.initialized)

    test_balance_equivalence(
        gpu.water, cpu.water; group = "water", rtol = 1.0f-3, atol = 5.0f-4,
    )
    test_balance_equivalence(
        gpu.nitrogen, cpu.nitrogen;
        group = "nitrogen", rtol = 2.0f-3, atol = 1.0f-3,
    )
    test_balance_equivalence(
        gpu.carbon, cpu.carbon;
        group = "carbon", rtol = 2.0f-3, atol = 1.0f-3,
        # The residual subtracts stocks of thousands of g C m-2. A 2 mg C m-2
        # CPU/GPU tolerance covers only final-bit Float32 cancellation.
        field_atol = (residual = 2.0f-3,),
    )
    test_balance_equivalence(
        gpu.thermal, cpu.thermal;
        group = "thermal", rtol = 2.0f-3, atol = 1.0f-1,
        # Energy residuals subtract O(1e8) J m-2 ledgers; keep state/flux fields
        # strict while allowing sub-joule cancellation in diagnostic closures.
        field_atol = (
            energy_residual = 1.0f0,
            percolation_energy_residual = 1.0f0,
        ),
    )
end
