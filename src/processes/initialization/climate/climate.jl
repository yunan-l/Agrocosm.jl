"""Daily weather forcing buffers on the active CPU/GPU backend."""
mutable struct DailyWeather{A}
    temp::A       # Daily mean air temperature (°C).
    prec::A       # Daily precipitation reaching the land surface (mm day⁻¹).
    swr::A        # Downward shortwave radiation (W m⁻²).
    lwr::A        # Net/downward longwave radiation forcing (W m⁻²).
    wind::A       # Near-surface wind speed (m s⁻¹).
    daily_co2::A  # Atmospheric CO₂ used on the current day (Pa).
    annual_co2::A # Annual atmospheric CO₂ forcing buffer (Pa).
end

"""Daily potential-evapotranspiration and radiation buffers."""
mutable struct PetPar{A}
    daylength::A # Astronomical daylight duration (h day⁻¹).
    par::A       # Daily incident photosynthetically active radiation (J m⁻² day⁻¹).
    eeq::A       # Equilibrium evapotranspiration demand (mm day⁻¹).
    albedo::A    # Effective land-surface albedo (fraction, 0–1).
end

"""Rolling climate buffers used by phenology and temperature-lag processes."""
mutable struct ClimBuf{A, M}
    temp::M       # Rolling daily air-temperature buffer (°C).
    mtemp::M      # Monthly mean air temperature (°C).
    mtemp20::M    # Long-term/20-year monthly mean temperature (°C).
    min_temp::M   # Recent minimum-temperature buffer used by phenology (°C).
    atemp::M      # Annual daily-temperature history (°C).
    atemp_mean::A # Annual mean air temperature (°C).
    V_req_a::A    # Acclimated vernalization requirement (day equivalent).
    V_req::A      # Current vernalization requirement (day equivalent).
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
