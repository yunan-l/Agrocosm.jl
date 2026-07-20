"""Daily weather forcing buffers on the active CPU/GPU backend."""
mutable struct DailyWeather{A}
    temp::A
    prec::A
    swr::A
    lwr::A
    wind::A
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

init_weather(cell_size::Int, device) = init_weather(Float32, cell_size, device)
function init_weather(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    cell_state() = device(zeros(T, cell_size))
    return DailyWeather(
        ntuple(_ -> cell_state(), 4)...,
        device(fill(T(lpjmlparams.volatil_wind), cell_size)),
        cell_state(),
        device(zeros(T, 1)),
    )
end

init_pet(cell_size::Int, device) = init_pet(Float32, cell_size, device)
function init_pet(::Type{T}, cell_size::Int, device) where {T <: AbstractFloat}
    cell_state() = device(zeros(T, cell_size))
    return PetPar(ntuple(_ -> cell_state(), 4)...)
end

init_climbuf(cell_size::Int, device; kwargs...) =
    init_climbuf(Float32, cell_size, device; kwargs...)
function init_climbuf(::Type{T},
                      cell_size::Int,
                      device;
                      NDAYS::Int = 31,
                      NMONTH::Int = 12,
                      NDAYS_YEAR::Int = 365,
                      n::Int = 5) where {T <: AbstractFloat}
    return ClimBuf(
        device(zeros(T, NDAYS, cell_size)),
        device(zeros(T, NMONTH, cell_size)),
        device(fill(T(-9999), NMONTH, cell_size)),
        device(zeros(T, n, cell_size)),
        device(zeros(T, NDAYS_YEAR, cell_size)),
        device(zeros(T, cell_size)),
        device(zeros(T, cell_size)),
        device(fill(T(-9999), cell_size)),
    )
end
