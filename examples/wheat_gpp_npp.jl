# # A ten-year wheat GPP/NPP simulation driven by real climate data
#
# This example reproduces the spirit of Agrocosm's original wheat example on the Terrarium stack: it
# runs a single-column temperate-cereals (wheat, CFT 1) crop for ten years, **driven by ten years of
# daily climate forcing** loaded from `climate_2000_2009.jld2`, sowing and harvesting each year through
# the crop management events, and recording the annual gross and net primary production, the peak leaf
# area index, and the harvest yield.
#
# The climate forcing is fed in the Terrarium-idiomatic way: `surface_climate_inputs` packs the daily
# series into Oceananigans `FieldTimeSeries` wrapped in Terrarium `FieldTimeSeriesInputSource`s, and
# `run!` interpolates them to the current time on every step — no manual per-step `set!` of the forcing.
#
# Run:  julia --project=. examples/wheat_gpp_npp.jl

using Agrocosm
using Terrarium
using JLD2

const NF = Float64
const SUBSTEPS = 144        # Δt = 600 s → 144 sub-steps per day
const Δt = 600.0

# --- load ten years of daily forcing for one grid cell -------------------------------------
climate = load(joinpath(@__DIR__, "climate_2000_2009.jld2"), "climate")
cell = 1
air_temperature = NF.(climate.temp[:, cell])         # °C
surface_shortwave = NF.(climate.swdown[:, cell])     # W/m²
ndays = length(air_temperature)                      # 3650 = 10 × 365
nyears = ndays ÷ 365

# One extra endpoint so the final day stays within the time-series range (no extrapolation).
times = (0:ndays) .* Terrarium.seconds_per_day(NF)
append_last(v) = vcat(v, v[end])

grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))
model = CropModel(grid, crop_pft("temperate cereals"))

# Time-varying temperature and shortwave come from the input sources; the remaining atmospheric forcing
# is held at physical constants (kept well-posed for the surface energy balance).
climate_inputs = surface_climate_inputs(
    grid, times;
    air_temperature = append_last(air_temperature),
    surface_shortwave_down = append_last(surface_shortwave),
)

integrator = initialize(model; inputs = climate_inputs, initializers = (temperature = 2.0,))
state = integrator.state
set!(state.air_pressure, 101325.0)
set!(state.specific_humidity, 0.006)
set!(state.windspeed, 2.0)
set!(state.surface_longwave_down, 300.0)
set!(state.CO2, 400.0)

calendar = CropCalendar(NF; sowing_day = 90, harvest_day = 240, residue_fraction = 0.25)

annual_gpp = zeros(nyears)
annual_npp = zeros(nyears)
annual_yield = zeros(nyears)
peak_lai = zeros(nyears)

# Step day by day: apply sowing/harvest at the season boundaries and advance a day; `run!` pulls the
# interpolated temperature/shortwave from the input sources automatically.
DAY = Terrarium.seconds_per_day(NF)
for year in 1:nyears
    for doy in 1:365
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

println("mean GPP (all years) = ", round(sum(annual_gpp) / nyears, digits = 4), " kgC/m²/yr")
