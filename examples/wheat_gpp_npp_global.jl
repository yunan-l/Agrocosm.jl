# # A global (multi-column) wheat GPP/NPP simulation
#
# The original Agrocosm example noted that the same process code runs over a *batch of independent grid
# cells* — one cell for a site run, a longer index vector for a larger batch. On the Terrarium stack
# that batch is a `ColumnRingGrid`: a soil/vegetation column at every cell of a global RingGrids grid,
# all advanced together. This example runs the wheat crop over such a grid, each column driven by its
# own daily climate series (the ten real climate cells tiled across the grid) and initialised with that
# cell's soil carbon from the initial-condition file, and maps the annual gross primary production.
#
# Global runs are heavy on CPU (cost scales with the number of columns × the fine surface timestep); the
# grid resolution and run length below are kept modest so the example finishes in a few minutes. The
# same script scales to finer grids by changing the grid resolution.
#
# Run:  julia --project=. examples/wheat_gpp_npp_global.jl

using Agrocosm
using Terrarium
using JLD2
import RingGrids
import Statistics: mean

const NF = Float64
const NYEARS = 1
const SUBSTEPS = 144
const Δt = 600.0

# A coarse global grid: a column at every cell of a full Gaussian grid (all cells active).
rings = RingGrids.FullGaussianGrid(2)
grid = ColumnRingGrid(CPU(), NF, ExponentialSpacing(Δz_max = 1.0, N = 20), rings)
lon, lat = RingGrids.get_lonlats(grid.rings)
ncolumns = length(lon)
println("global grid: ", ncolumns, " columns")

# Initial-condition file (per cell): management + initial soil carbon. Tile the available cells across
# the grid columns.
ic = load_crop_initial_conditions(joinpath(@__DIR__, "initial_wheat.jld2"))
ncells = length(ic.sowing_day)
cell_of(c) = mod1(c, ncells)
pft = crop_pft("temperate cereals")

# Load the daily forcing and tile the climate cells across the grid columns (a `time × column` matrix).
climate = load(joinpath(@__DIR__, "climate_2000_2009.jld2"), "climate")
ndays = NYEARS * 365
tile(v) = NF[v[t, cell_of(c)] for t in 1:ndays, c in 1:ncolumns]
with_endpoint(m) = vcat(m, m[end:end, :])
temp = with_endpoint(tile(climate.temp))
shortwave = with_endpoint(tile(climate.swdown))
times = (0:ndays) .* Terrarium.seconds_per_day(NF)

# The phenology/heat-unit requirement and management calendar are scalar (per-model), so use a
# representative value from the initial conditions; the initial soil carbon is set per column below.
model = CropModel(
    grid, pft;
    vegetation = CropVegetation(NF, pft; heat_unit_requirement = mean(ic.heat_unit_requirement)),
)
climate_inputs = surface_climate_inputs(
    grid, times; air_temperature = temp, surface_shortwave_down = shortwave,
)

integrator = initialize(model; inputs = climate_inputs, initializers = (temperature = 2.0,))
state = integrator.state
set!(state.air_pressure, 101325.0)
set!(state.specific_humidity, 0.006)
set!(state.windspeed, 2.0)
set!(state.surface_longwave_down, 300.0)
set!(state.CO2, 400.0)

# Per-column initial soil carbon from the (tiled) initial-condition cells.
fast_carbon = interior(state.fast_carbon)
slow_carbon = interior(state.slow_carbon)
for c in 1:ncolumns
    fast_carbon[c, 1, :] .= ic.fast_carbon[cell_of(c)]
    slow_carbon[c, 1, :] .= ic.slow_carbon[cell_of(c)]
end

sowing_day = ic.sowing_day[1]
calendar = CropCalendar(NF; sowing_day, harvest_day = mod1(sowing_day + 250, 365),
    residue_fraction = mean(ic.residue_fraction))

DAY = Terrarium.seconds_per_day(NF)
annual_gpp = zeros(ncolumns)
peak_lai = zeros(ncolumns)

println("running ", NYEARS, " year(s) over ", ncolumns, " columns...")
for year in 1:NYEARS
    for doy in 1:365
        doy == calendar.sowing_day && sow!(integrator, calendar)
        doy == calendar.harvest_day && harvest!(integrator, calendar)
        run!(integrator; steps = SUBSTEPS, Δt = Δt)
        gpp = interior(state.gross_primary_production)[:, 1, 1]
        lai = interior(state.leaf_area_index)[:, 1, 1]
        annual_gpp .+= max.(0.0, gpp) .* DAY
        peak_lai .= max.(peak_lai, lai)
    end
    println("  year ", year, " done")
    flush(stdout)
end

println("per-column annual GPP (kgC/m²/yr): min=", round(minimum(annual_gpp), digits = 4),
    " mean=", round(mean(annual_gpp), digits = 4), " max=", round(maximum(annual_gpp), digits = 4))
println("per-column peak LAI:              min=", round(minimum(peak_lai), digits = 2),
    " mean=", round(mean(peak_lai), digits = 2), " max=", round(maximum(peak_lai), digits = 2))
println("global-mean annual GPP = ", round(mean(annual_gpp), digits = 4), " kgC/m²/yr")
