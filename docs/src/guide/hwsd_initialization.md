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
`soil_initial_state` creates the seven required initialization arrays from
the source-neutral targets:

```julia
grid = read_grid("grid.nc")
soil = read_soil_data(catalog, grid; selection)
targets = read_soil_cn_targets("hwsd_cn_targets.nc")
initial_state = soil_initial_state(targets, soil)
```

The initialization defaults are:

- `fastc = 0.4 × SOC` and `slowc = 0.6 × SOC`;
- organic N uses the same 40:60 split after reserving the initial mineral N;
- surface, incorporated, and root litter C/N start at zero;
- soil liquid water starts at field capacity using Agrocosm's pedotransfer
  equations and the HWSD SOC stock.

For layer ``l``, let ``f=0.4``, ``s=0.6``, and let ``m=0.01`` be the LPJmL
fresh-soil mineral-N fraction. The initialized pools are

```math
C_{fast,l}=fC_{SOC,l},\qquad C_{slow,l}=sC_{SOC,l}.
```

Each of the nitrate and ammonium pools is initialized to ``m`` times slow
organic N. To conserve the HWSD total-N target exactly, Agrocosm first computes

```math
N_{org,l}=\frac{N_{HWSD,l}}{1+2ms},
```

and then assigns

```math
N_{fast,l}=fN_{org,l},\qquad
N_{slow,l}=sN_{org,l},\qquad
N_{NO_3,l}=N_{NH_4,l}=mN_{slow,l}.
```

Consequently,

```math
N_{fast,l}+N_{slow,l}+N_{NO_3,l}+N_{NH_4,l}=N_{HWSD,l}.
```

The 40:60 split matches the mean partition of the legacy ten-cell spun-up
fixture. It is only a reproducible starting guess, not an HWSD observation or
an equilibrium claim.

## Target-constrained agricultural warm-up

Production initialization does not use the raw 40:60 state directly. Agrocosm
runs an adaptive agricultural warm-up before the first reported day. Ten years
is the minimum duration; cells that have not stabilized cause the complete
batch to continue, up to a configured maximum.

The initial HWSD stocks define fixed layer targets:

```math
C^*_{l}=C_{fast,l}^{(0)}+C_{slow,l}^{(0)}=C_{SOC,l},
```

```math
N^*_{l}=N_{fast,l}^{(0)}+N_{slow,l}^{(0)}+
N_{NO_3,l}^{(0)}+N_{NH_4,l}^{(0)}=N_{HWSD,l}.
```

After each warm-up year, fast and slow C are scaled proportionally so their
sum returns to ``C^*_l``. Fast and slow organic N plus nitrate and ammonium are
scaled proportionally so their sum returns to ``N^*_l``. If a current layer
sum is numerically zero, the original initialized partition is restored.

Fresh litter is deliberately excluded from these constraints because HWSD SOC
and total N describe mineral-soil stocks rather than an above-soil litter
layer. Litter therefore develops from crop residues, roots, decomposition,
tillage, and management without being reset. The annual target correction

```math
\Delta C_{target,l}=C^*_l-C_{l}^{pre},\qquad
\Delta N_{target,l}=N^*_l-N_{l}^{pre}
```

is an initialization diagnostic, not a production carbon or nitrogen flux.
It measures how strongly the unconstrained annual processes would move the
HWSD target pools.

For each cell, convergence requires all of the following for a configured
number of consecutive years:

- annual relative change in total soil C and N, including litter, is below the
  relative tolerance;
- annual change in fast-C and fast-N fractions is below the pool-fraction
  tolerance;
- the absolute annual target correction, relative to the initial cell total,
  is below the relative tolerance.

The global CPU workflow cycles 1901--1930 climate five times and runs every
selected cell for 150 years. It still reports whether each cell passes for
three consecutive years. Keeping all cells on the same fixed calendar makes
the checkpoint deterministic and exposes the residual unconverged fraction.

Recommended procedure:

1. Build the HWSD state and initialize water at field capacity.
2. Run the 150-year target-constrained agricultural warm-up, using
   the same crop, fertilization, irrigation, residue, and tillage configuration
   as the target experiment.
3. Use observed historical forcing when available. If forcing must be cycled,
   pass a complete multi-year block rather than one anomalous year.
4. Evaluate convergence after the five complete 30-year forcing cycles.
5. Discard warm-up outputs but save the final native Agrocosm checkpoint.
6. Use a separate diagnostic run when daily C, N, water, and energy closure
   must be audited; the production warm-up deliberately does not allocate
   daily balance ledgers.
7. Inspect annual target corrections, litter, fast/slow fractions, mineral N,
   total C/N, the converged-cell fraction, and unresolved cells.

Reaching the maximum without convergence does not silently accept the state.
The report uses `target_constrained_maximum_years`, and production output from
that checkpoint is treated as a baseline pending review of fallback soils,
forcing, management, thresholds, or genuinely slow pool dynamics.

## Provenance and reproducibility

Every generated file records the HWSD version, preprocessing version, vertical
rule, component policy, spatial policy, coverage, uncertainty, and donor
metadata when applicable. Production runs should retain the preprocessing
configuration and warm-up checkpoint alongside model outputs.

## Running the warm-up

Initialize the production simulation normally, then warm its state before the
first reported day:

```julia
simulation = initialize_simulation(
    cft1,
    initial_data;
    days = size(production_climate.temp, 1),
    diagnostics = false,
)

report = agricultural_warmup!(
    simulation,
    warmup_climate;
    years = 10,
    maximum_years = 100,
    target_constrained = true,
    consecutive_years = 3,
    relative_tolerance = 0.01,
    pool_fraction_tolerance = 0.01,
    required_converged_fraction = 1.0,
)
@assert simulation.simulated_days == 0
run_simulation!(simulation, production_climate; spinup = false)
```

`warmup_climate` must contain one or more complete 365-day years. If fewer
years are supplied than the warm-up needs, the years cycle in their original
order. The report stores the actual duration, convergence status, annual target
corrections, converged-cell fractions, and host-side matrices for total soil
C/N, litter, fast and slow C/N, mineral N, and soil water. Warm-up outputs and
balance ledgers are not mixed with the production run.

The server production script exposes the same controls:

```toml
[run]
warmup_target_constrained = true
warmup_minimum_years = 150
warmup_maximum_years = 150
warmup_consecutive_years = 3
warmup_relative_tolerance = 0.01
warmup_pool_fraction_tolerance = 0.01
warmup_required_converged_fraction = 1.0
```

The global production configuration cycles the 1901--1930 climate five times,
for a fixed 150-year target-constrained warm-up. Convergence compares states at
the same position in the 30-year forcing cycle. HWSD mineral-soil C and total-N
targets remain fixed while litter and the fast/slow pool allocation develop
under the model processes.

For the complete canonical grid, run the bounded-memory HWSD raster pipeline
on the server:

```bash
julia --project=lib/AgrocosmData \
  lib/AgrocosmData/scripts/prepare_canonical_hwsd.jl \
  /absolute/path/HWSD2 \
  /absolute/path/grid.nc \
  /absolute/path/hwsd_canonical.nc \
  /absolute/path/hwsd_canonical_qc.toml
```

The output contains every valid canonical `cellid`, five-layer SOC and total
N, coverage and uncertainty, fallback donor identifiers and distances, and
soil area. The QC report records unresolved cells and layer-wise C/N
conservation errors. The command fails if fallback cannot resolve every cell
or conservation exceeds numerical tolerance. `mdb-sql` from MDB Tools must be
available because HWSD v2.01 distributes its mapping-unit table as an MDB file.
