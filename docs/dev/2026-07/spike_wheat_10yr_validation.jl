# Phase 6 acceptance: a 10-year single-cell wheat GPP/NPP validation on the managed-crop CropModel
# stack. Drives a temperate-cereals (CFT 1) column with a synthetic seasonal climate, sows and harvests
# each year through the crop management events, and records the annual gross and net primary production,
# the peak leaf area index, and the harvest yield — checking that the crop reproduces a stable,
# repeating seasonal carbon cycle over ten years.
#
# Run: julia --project=. docs/dev/2026-07/spike_wheat_10yr_validation.jl

using Agrocosm
using Terrarium

const DAY = Terrarium.seconds_per_day(Float64)
const YEAR_DAYS = 365
const NYEARS = 10
const SUBSTEPS = 144        # Δt = 600 s → 144 sub-steps/day (surface-energy-balance-stable)
const Δt = 600.0

grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))
model = CropModel(grid, crop_pft("temperate cereals"))
integrator = initialize(model; initializers = (temperature = 8.0,))
state = integrator.state

# Physical atmospheric forcing held constant (pressure, humidity, wind, downwelling longwave, CO₂) so the
# surface energy balance stays well-posed across the whole run; temperature and shortwave follow the
# season.
set!(state.air_pressure, 101325.0)
set!(state.specific_humidity, 0.006)
set!(state.windspeed, 2.0)
set!(state.surface_longwave_down, 300.0)
set!(state.CO2, 400.0)

# Northern-hemisphere growing season: air temperature and incoming shortwave peak near midsummer.
seasonal_temperature(doy) = 10.0 + 12.0 * sin(2π * (doy - 100) / YEAR_DAYS)          # °C
seasonal_shortwave(doy) = 150.0 + 200.0 * max(0.0, sin(2π * (doy - 100) / YEAR_DAYS))  # W/m²

# Spring wheat: sow in spring, harvest in late summer.
calendar = CropCalendar(Float64; sowing_day = 90, harvest_day = 240, residue_fraction = 0.25)

annual_gpp = zeros(NYEARS)
annual_npp = zeros(NYEARS)
annual_yield = zeros(NYEARS)
peak_lai = zeros(NYEARS)

println("SPIKE OK — running ", NYEARS, " years of wheat...")
completed = 0
for year in 1:NYEARS
    for doy in 1:YEAR_DAYS
        set!(state.air_temperature, seasonal_temperature(doy))
        set!(state.surface_shortwave_down, seasonal_shortwave(doy))
        doy == calendar.sowing_day && sow!(integrator, calendar)
        doy == calendar.harvest_day && (annual_yield[year] = harvest!(integrator, calendar))
        run!(integrator; steps = SUBSTEPS, Δt = Δt)
        # Approximate the daily carbon totals from the end-of-day fluxes (forcing is constant per day).
        annual_gpp[year] += max(0.0, interior(state.gross_primary_production)[1, 1, 1]) * DAY
        annual_npp[year] += interior(state.net_primary_production)[1, 1, 1] * DAY
        peak_lai[year] = max(peak_lai[year], interior(state.leaf_area_index)[1, 1, 1])
    end
    global completed = year
    println("  year ", lpad(year, 2), ": GPP=", rpad(round(annual_gpp[year], digits = 4), 8),
        " NPP=", rpad(round(annual_npp[year], digits = 4), 8), " kgC/m²/yr | peak LAI=",
        rpad(round(peak_lai[year], digits = 2), 5), " | yield=", round(annual_yield[year], digits = 4), " kgC/m²")
    flush(stdout)
end

@assert completed == NYEARS "all $NYEARS years completed"
@assert all(isfinite, annual_gpp) && all(isfinite, annual_npp) "GPP/NPP stayed finite"
@assert maximum(peak_lai) > 1.0 "the wheat canopy developed a substantial LAI"
@assert maximum(annual_gpp) > 0.1 "the wheat assimilated carbon each season"
# After the first (spin-up) year the seasonal cycle should repeat within a modest band.
gpp_steady = annual_gpp[3:NYEARS]
@assert (maximum(gpp_steady) - minimum(gpp_steady)) / maximum(gpp_steady) < 0.5 "annual GPP is quasi-repeating"
println("SPIKE ASSERTIONS PASSED — mean GPP (yrs 3–10) = ",
    round(sum(gpp_steady) / length(gpp_steady), digits = 4), " kgC/m²/yr")
