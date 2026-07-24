# # A ten-year wheat GPP/NPP simulation
#
# This example runs a single-column temperate-cereals (wheat, CFT 1) crop for ten years under a
# synthetic seasonal climate, sowing and harvesting each year through the crop management events, and
# records the annual gross and net primary production, the peak leaf area index, and the harvest yield.
# It exercises the whole managed-crop stack — phenology, photosynthesis, the carbon and nitrogen pools,
# the soil carbon–nitrogen biogeochemistry, and the sowing/harvest lifecycle — and shows that the crop
# settles into a stable, repeating seasonal carbon cycle.
#
# Run:  julia --project=. examples/wheat_gpp_npp.jl

using Agrocosm
using Terrarium

const DAY = Terrarium.seconds_per_day(Float64)
const YEAR_DAYS = 365
const NYEARS = 10
const SUBSTEPS = 144        # Δt = 600 s → 144 sub-steps per day
const Δt = 600.0

# A single wheat column.
grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))
model = CropModel(grid, crop_pft("temperate cereals"))
integrator = initialize(model; initializers = (temperature = 8.0,))
state = integrator.state

# Physical atmospheric forcing held constant (pressure, humidity, wind, downwelling longwave, CO₂) so
# the surface energy balance stays well-posed; temperature and shortwave follow the season.
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

# Step day by day: update the daily forcing, apply sowing/harvest at the season boundaries, advance the
# continuous dynamics for a day, and accumulate the daily carbon fluxes.
for year in 1:NYEARS
    for doy in 1:YEAR_DAYS
        set!(state.air_temperature, seasonal_temperature(doy))
        set!(state.surface_shortwave_down, seasonal_shortwave(doy))
        doy == calendar.sowing_day && sow!(integrator, calendar)
        doy == calendar.harvest_day && (annual_yield[year] = harvest!(integrator, calendar))
        run!(integrator; steps = SUBSTEPS, Δt = Δt)
        annual_gpp[year] += max(0.0, interior(state.gross_primary_production)[1, 1, 1]) * DAY
        annual_npp[year] += interior(state.net_primary_production)[1, 1, 1] * DAY
        peak_lai[year] = max(peak_lai[year], interior(state.leaf_area_index)[1, 1, 1])
    end
    println("year ", lpad(year, 2), ": GPP=", rpad(round(annual_gpp[year], digits = 4), 8),
        " NPP=", rpad(round(annual_npp[year], digits = 4), 8), " kgC/m²/yr | peak LAI=",
        rpad(round(peak_lai[year], digits = 2), 5), " | yield=", round(annual_yield[year], digits = 4), " kgC/m²")
    flush(stdout)
end

# After the first (spin-up) year the crop reaches a stable, repeating seasonal cycle.
mean_gpp = sum(annual_gpp[3:NYEARS]) / (NYEARS - 2)
println("mean GPP (years 3–10) = ", round(mean_gpp, digits = 4), " kgC/m²/yr")
