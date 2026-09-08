# Preparing global CFT input subsets

The bounded preprocessing script extracts the 12 canonical crops in both rainfed and
irrigated form without changing the spatial grid. A management-only
configuration can preserve every available source year for both transient
and fixed-year management. Climate and CO₂ extraction
remain optional for compact test fixtures.

Land use, fertilizer, and manure use 24 selected source bands. Sowing date and
PHU use their 24 rainfed/irrigated bands directly. Residue management provides
12 crop bands shared by rainfed and irrigated patches. The exact source band
positions are declared explicitly through `cft_indices` in the configuration.

## 1. Configure input and output paths

Copy the example configuration:

```bash
cp lib/AgrocosmData/config/global_cft_subset.example.toml \
   global_cft_subset.toml
```

Edit the input paths and `output_directory`. Each entry explicitly declares
the NetCDF variable because filenames and variable names are independent.

| Dataset | Selected output bands | Default variable |
|---|---:|---|
| land use | 24 | `landfrac` |
| fertilizer | 24 | `fertilizer` |
| manure | 24 | `manure` |
| residue on field | 12 | `residuefrac` |
| sowing date | 24 | `sdate` |
| PHU | 24 | `phusum` |

The CFT dimension may be named `cft`, `cft`, `crop`, or `band`, and it may
occur in any dimension position. The script requires exactly one such
dimension in each management variable.

## 2. Run the extraction

From the Agrocosm repository root:

```bash
julia --project=lib/AgrocosmData \
  -e 'import Pkg; Pkg.instantiate()'

julia --project=lib/AgrocosmData \
  lib/AgrocosmData/scripts/prepare_global_cft_subset.jl \
  /absolute/path/global_cft_subset.toml
```

The script performs bounded time-chunk reads. `chunk_length = 31` means that at
most 31 records along the time dimension are copied in one I/O batch. It does
not change values or the selected time range. Existing output files are not
overwritten.

Set `management_years = "all"` to preserve every available management year.
At runtime, `management.mode = "transient"` selects the corresponding year,
while `management.mode = "fixed"` repeats the configured `fixed_year`.
Static inputs such as the current sowing-date file remain static.

When a `[climate]` section is present, every selected year must contain exactly
365 rows. A file containing 29 February or any other annual row count is
rejected because Agrocosm uses a no-leap 365-day calendar.

## 3. Resulting files

With the supplied management-only template, the output directory contains six
management NetCDF files. NetCDF dimension order, coordinate variables, data
type, fill value, units, calendar, and source attributes are preserved.
Additional global attributes record:

- the absolute source path;
- the selected variable;
- the selected source CFT indices;
- selected management years when a time subset is requested.

The CFT dimension remains present with 24 entries for rainfed/irrigated
management and 12 entries for shared residue management. Output coordinates
are renumbered consecutively; the source positions remain recorded in global
attributes.

## 4. Validate before transfer

Check the generated file headers and sizes:

```bash
ncdump -h /output/path/landuse.nc
ncdump -h /output/path/phu.nc
```

Expected properties:

- management CFT/band dimension: `24`, or `12` for residue management;
- longitude and latitude dimensions: unchanged;
- management time dimension: unchanged when `management_years = "all"`;
- source CFT indices recorded in `agrocosm_source_cft_indices`;
- all requested source years retained.

Then point the AgrocosmData catalog at the extracted files. The catalog maps
the 12 rainfed and 12 irrigated output positions through the canonical CFT
registry; residue positions are shared between both water-management modes.

## 5. Scope

This extraction reduces source-band volume, not space. The historical
management outputs support both fixed-year and transient-management
experiments on the real `720 × 280` grid. Any optional short climate subset is
only a portable test fixture, not a scientifically complete experiment.

## 6. CPU production workflow

Copy and edit the production configuration:

```bash
cp scripts/global_wheat_cpu.example.toml global_wheat_cpu.toml

julia --project=. scripts/run_global_wheat_cpu.jl \
  /absolute/path/global_wheat_cpu.toml
```

The `management.mode = "fixed"` option repeats the configured `fixed_year` throughout
the run. With `management.mode = "transient"`, each simulation year reads its
corresponding management row. Years before the file begins repeat its first
row, and years after it ends repeat its last row, matching LPJmL's boundary
behavior. Annual PHU is installed only when a new crop is sown, so a winter
crop already growing across 1 January retains the PHU assigned in its sowing
year. `landfrac > 0` selects crop cells but never scales single-cell processes.

The climate files must cover the configured warm-up and production periods.
Warm-up cycles `warmup_climate_start_year` through `warmup_climate_end_year`;
its duration and management policy are explicit run settings. Convergence
compares soil states at the same position in the forcing cycle, not unlike
adjacent climate years. Choose periods and acceptance criteria in the external
experiment configuration.

The runner performs the configured streamed agricultural warm-up, writes and
exactly restores a native warm-up checkpoint, checkpoints after the first
production year, restores, and continues through `simulation_end_year`.
Annual compact outputs are streamed to NetCDF and reconstructed onto the
canonical longitude/latitude grid.

The full domain runs with daily balance ledgers disabled. A configurable small
canonical subset repeats the same warm-up and production run with diagnostics to
produce sampled C/N/water/energy closure. The output directory also contains a
pre-allocation memory estimate and a warm-up C/N drift report.

## 7. CUDA production workflow

The CUDA entry uses the same backend-neutral AgrocosmData inputs and production
contract. HWSD C/N, crop management, soil properties, and each bounded climate
block are transferred during initialization or immediately before execution;
checkpoint and streamed output arrays are copied back to the host.

First validate the node and the backend-neutral initialization path:

```bash
julia --project=. -e 'using CUDA; @assert CUDA.functional(); CUDA.versioninfo()'
julia --project=. test/simulations/test_global_initialization_gpu.jl
julia --project=. test/simulations/test_daily_crop_C3_gpu.jl
```

Use a separate output directory from the CPU baseline, then run:

```bash
JULIA_NUM_THREADS=4 julia --project=. \
  scripts/run_global_wheat_gpu.jl \
  /absolute/path/global_wheat_gpu.toml
```

Set `run.device_id` to the intended visible CUDA device.
Review `recommended_device_peak_gib` in `memory_preflight.toml` before the full
run. Validate a bounded subset before setting `cell_limit = 0`.
