"""Daily weather forcing buffers on the active CPU/GPU backend."""
mutable struct DailyWeather{A}
    temp::A
    prec::A
    swr::A
    lwr::A
    daily_co2::A
    annual_co2::A
end

"""Daily potential-evapotranspiration and radiation buffers."""
mutable struct PetPar{A}
    daylength::A
    par::A
    eeq::A
    albedo::A
end

"""Rolling climate buffers used by phenology and temperature-lag processes."""
mutable struct ClimBuf{A, M}
    temp::M
    mtemp::M
    mtemp20::M
    min_temp::M
    atemp::M
    atemp_mean::A
    V_req_a::A
    V_req::A
end

function init_weather(cell_size::Int, device)
    cell_state() = device(zeros(Float32, cell_size))
    return DailyWeather(
        ntuple(_ -> cell_state(), 5)...,
        device(zeros(Float32, 1)),
    )
end

function init_pet(cell_size::Int, device)
    cell_state() = device(zeros(Float32, cell_size))
    return PetPar(ntuple(_ -> cell_state(), 4)...)
end

function init_climbuf(cell_size::Int,
                      device;
                      NDAYS::Int = 31,
                      NMONTH::Int = 12,
                      NDAYS_YEAR::Int = 365,
                      n::Int = 5)
    return ClimBuf(
        device(zeros(Float32, NDAYS, cell_size)),
        device(zeros(Float32, NMONTH, cell_size)),
        device(fill(-9999.0f0, NMONTH, cell_size)),
        device(zeros(Float32, n, cell_size)),
        device(zeros(Float32, NDAYS_YEAR, cell_size)),
        device(zeros(Float32, cell_size)),
        device(zeros(Float32, cell_size)),
        device(fill(-9999.0f0, cell_size)),
    )
end
