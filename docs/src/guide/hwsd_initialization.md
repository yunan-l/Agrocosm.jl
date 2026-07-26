# HWSD soil initialization

Agrocosm can initialize soil water, carbon, and nitrogen without an LPJmL
restart. The workflow uses HWSD v2.01 for soil organic carbon (SOC) and total
nitrogen targets, then constructs the model pools before the production run.
HWSD preprocessing is an offline CPU task and is not part of the differentiable
daily transition.

## Spatial alignment

`grid.nc` is the canonical grid for climate, management, soil, and HWSD data.
Every lookup uses the one-based `(longitude_index, latitude_index)` stored by
Agrocosm; code must not reconstruct coordinates from an assumed origin. The
current grid has 720 longitudes and 280 latitudes at 0.5° resolution.

HWSD is sampled from its 30 arc-second raster. The 3600 source pixels inside
each 0.5° target cell are aggregated by spherical pixel area. Raster NODATA
and mapping units explicitly classified as water (`WR`) or glacier (`GG`) are
excluded from the soil-average denominator rather than counted as missing
soil. This makes coastal-cell values representative of their soil-covered
portion.

## Concentrations to stocks

HWSD provides seven layers: 0–20, 20–40, 40–60, 60–80, 80–100, 100–150,
and 150–200 cm. For each soil component and layer, Agrocosm converts organic
carbon, total nitrogen, bulk density, and coarse fragments to area stocks:

```text
SOC [gC m⁻²] = OC [wt%] × bulk density [g cm⁻³]
                × thickness [cm] × fine-earth fraction × 100

N [gN m⁻²] = total N [g kg⁻¹] × bulk density [g cm⁻³]
              × thickness [cm] × fine-earth fraction × 10
```

The fine-earth fraction is `1 - coarse_fragments / 100`.

## Soil components and shallow profiles

An HWSD mapping unit can contain several soil components with different
`SHARE` and `ROOT_DEPTH` values. The root-depth classes are:

| Code | Meaning | Last represented HWSD layer |
|---:|---|---:|
| 1 | deep, greater than 100 cm | D7 |
| 2 | moderately deep, less than 100 cm | D5 |
| 3 | shallow, less than 50 cm | D2 |
| 4 | very shallow, less than 10 cm | D1 |

Missing attributes below the last represented layer are structural zeros: the
component contributes no C or N stock there. Its original `SHARE` remains in
the weighted mean. Available deep components are never renormalized to 100%,
which would overestimate deep stocks in mixed shallow/deep mapping units.

A missing value inside the represented root depth is a true data gap. A layer
is accepted only when at least 99% of component share is resolved; otherwise
the target remains missing and is eligible for an explicitly recorded donor
fallback. Structural-zero layers are marked uncertain because root depth is a
categorical constraint rather than an exact measured boundary.

## Vertical mapping

Stocks are mapped by exact depth overlap to Agrocosm's layers:

```text
0–0.2, 0.2–0.5, 0.5–1.0, 1.0–2.0, and 2.0–3.0 m
```

HWSD ends at 2 m. By default, the 1.5–2 m stock density is extended through
2–3 m and that layer is marked uncertain. Use `deep_rule=:missing` when this
extrapolation is unsuitable.

## Missing-profile fallback

Nearest-cell filling is a last resort, not the treatment for shallow soil or
coastal pixels. It is used only when a true HWSD gap leaves an unresolved model
layer. Filled profiles record the donor longitude, latitude, distance, original
minimum coverage, and uncertainty in the output NetCDF. Preprocessing fails if
no complete donor exists inside the configured search radius.

## Constructing model state

The preprocessing output contains SOC and total-N targets, not model pools.
`hwsd_initial_state` creates the seven required initialization arrays:

```julia
grid = read_grid("grid.nc")
soil = read_soil_data(catalog, grid; selection)
targets = read_soil_cn_targets("hwsd_cn_targets.nc")
initial_state = hwsd_initial_state(targets, soil)
```

The current no-spin-up defaults are:

- `fastc = 0.4 × SOC` and `slowc = 0.6 × SOC`;
- organic N uses the same 40:60 split;
- organic N is reduced slightly so organic plus the default initial nitrate and
  ammonium pools conserve the HWSD total-N target;
- surface, incorporated, and root litter C/N start at zero;
- soil liquid water starts at field capacity using Agrocosm's pedotransfer
  equations and the HWSD SOC stock.

The 40:60 split matches the mean partition of the legacy ten-cell spun-up
fixture. It is an initialization assumption, not an HWSD observation.

## Recommended ten-year warm-up

A ten-year pre-run is reasonable for generating crop residue and allowing
litter and faster soil pools to adapt to local climate and management. It
should be called a **warm-up**, not a complete soil-carbon spin-up: ten years is
too short for the slow SOC pool to reach equilibrium.

Recommended procedure:

1. Build the HWSD state and initialize water at field capacity.
2. Run ten years before the reported simulation, using the same crop,
   fertilization, irrigation, residue, and tillage configuration as the target
   experiment.
3. Use observed historical forcing when available. If forcing must be cycled,
   repeat a multi-year block rather than one anomalous year.
4. Discard warm-up outputs but save the final native Agrocosm checkpoint.
5. Verify daily C, N, water, and energy closure throughout the warm-up.
6. Report annual litter C/N, fast C/N, mineral N, heterotrophic respiration,
   NPP, and yield. The last three-year mean should not show strong monotonic
   drift relative to the preceding three years.

Ten years is accepted when litter and fast pools stabilize and crop behaviour
is plausible. Failure to stabilize means the warm-up must be lengthened or the
initial partition reconsidered. Slow-pool convergence is deliberately not an
acceptance criterion for this interim workflow; a future full spin-up may need
many decades or longer and should be treated as a separate offline stage.

## Provenance and reproducibility

Every generated file records the HWSD version, preprocessing version, vertical
rule, component policy, spatial policy, coverage, uncertainty, and donor
metadata when applicable. Production runs should retain the preprocessing
configuration and warm-up checkpoint alongside model outputs.
