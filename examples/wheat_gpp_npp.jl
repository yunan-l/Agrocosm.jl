# # A ten-year wheat GPP/NPP simulation driven by real climate + initial-condition data
#
# This reproduces the spirit of Agrocosm's original wheat example on the Terrarium stack. A single
# temperate-cereals (wheat, CFT 1) column is run for ten years, **driven by ten years of daily climate
# forcing** and **initialised from the site's initial-condition file** — the sowing date, the
# phenological-heat-unit requirement, the residue fraction, and the initial soil carbon pools all come
# from `initial_wheat.jld2`. It records the annual gross and net primary production, the peak leaf area
# index, and the harvest yield.
#
# The forcing is fed the Terrarium-idiomatic way (`surface_climate_inputs` → `FieldTimeSeries` input
# sources, interpolated each step); the site setup is read with `load_crop_initial_conditions`.
#
# Run:  julia --project=. examples/wheat_gpp_npp.jl

using Agrocosm
using Terrarium
using JLD2

const NF = Float64
const SUBSTEPS = 144        # Δt = 600 s → 144 sub-steps per day
const Δt = 600.0
const cell = 1              # which input grid cell to run

# --- site setup from the initial-condition file --------------------------------------------
ic = load_crop_initial_conditions(joinpath(@__DIR__, "initial_wheat.jld2"))
sowing_day = ic.sowing_day[cell]
harvest_day = mod1(sowing_day + 250, 365)   # ~250-day season (autumn-sown winter wheat here)
pft = crop_pft("temperate cereals")

# --- ten years of daily forcing for the same cell ------------------------------------------
climate = load(joinpath(@__DIR__, "climate_2000_2009.jld2"), "climate")
air_temperature = NF.(climate.temp[:, cell])         # °C
surface_shortwave = NF.(climate.swdown[:, cell])     # W/m²
ndays = length(air_temperature)                      # 3650 = 10 × 365
nyears = ndays ÷ 365
times = (0:ndays) .* Terrarium.seconds_per_day(NF)   # one extra endpoint (no extrapolation)
append_last(v) = vcat(v, v[end])

grid = ColumnGrid(CPU(), ExponentialSpacing(Δz_max = 1.0, N = 20))

# Model configured from the initial conditions: the site's heat-unit requirement (phenology) and its
# initial fast/slow soil carbon pools.
model = CropModel(
    grid, pft;
    vegetation = CropVegetation(NF, pft; heat_unit_requirement = ic.heat_unit_requirement[cell]),
    soil_biogeochemistry = CropSoilBiogeochemistry(NF;
        initial_fast_carbon = ic.fast_carbon[cell], initial_slow_carbon = ic.slow_carbon[cell]),
)

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

calendar = CropCalendar(NF; sowing_day, harvest_day, residue_fraction = ic.residue_fraction[cell])
println("site: sowing day ", sowing_day, ", harvest day ", harvest_day,
    ", PHU ", ic.heat_unit_requirement[cell], " °C·days, residue ", round(ic.residue_fraction[cell], digits = 2))

annual_gpp = zeros(nyears)
annual_npp = zeros(nyears)
annual_yield = zeros(nyears)
peak_lai = zeros(nyears)

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
