# Preparing the global wheat test subset

The production source files are intentionally kept on the server. A small
preprocessing script creates a portable global rainfed-wheat dataset without
changing the spatial grid:

- management files retain every longitude, latitude, and time record but only
  the first rainfed crop band;
- climate files retain every grid cell but only the first two complete calendar
  years shared by temperature, precipitation, longwave, and shortwave;
- the annual CO₂ text file is reduced to the same two years.

The first LPJmL crop position is **temperate cereals**. Agrocosm uses its wheat
parameter set for the initial production test. In all current 64/32/24/16-band
management products, its rainfed position is band 1. The script ignores the
irrigated bands.

## 1. Configure server paths

Copy the example configuration:

```bash
cp lib/AgrocosmData/config/global_wheat_subset.example.toml \
   global_wheat_subset.toml
```

Edit the input paths and `output_directory`. Each entry explicitly declares
the NetCDF variable because filenames and variable names are independent.

| Dataset | Source bands | Rainfed wheat band | Default variable |
|---|---:|---:|---|
| land use | 64 | 1 | `landfrac` |
| fertilizer | 32 | 1 | `fertilizer` |
| manure | 32 | 1 | `manure` |
| residue on field | 16 | 1 | `residuefrac` |
| sowing date | 24 | 1 | `sdate` |
| PHU | 24 | 1 | `phusum` |

The CFT dimension may be named `cft`, `cft`, `crop`, or `band`, and it may
occur in any dimension position. The script requires exactly one such
dimension in each management variable.

## 2. Run the extraction

From the Agrocosm repository root:

```bash
julia --project=lib/AgrocosmData \
  -e 'import Pkg; Pkg.instantiate()'

julia --project=lib/AgrocosmData \
  lib/AgrocosmData/scripts/prepare_global_wheat_subset.jl \
  /absolute/path/global_wheat_subset.toml
```

The script performs bounded time-chunk reads. `chunk_length = 31` means that a
large daily field is copied at most 31 source days at a time; it never loads
the complete climate file. Existing output files are not overwritten.

For climate, the script decodes the CF time coordinate, finds the first two
complete January–December years in the reference file, and requires those same
years in every other climate file. Every selected year must contain exactly
365 rows. A file containing 29 February or any other annual row count is
rejected because the current production forcing already uses a no-leap
365-day calendar.

## 3. Resulting files

The output directory contains six management NetCDF files, four climate
NetCDF files, and a two-row CO₂ text file. NetCDF dimension order, coordinate
variables, data type, fill value, units, calendar, and source attributes are
preserved. Additional global attributes record:

- the absolute source path;
- the selected variable;
- rainfed CFT index 1 for management data;
- selected climate years for climate data.

The CFT dimension remains present with length one. Keeping it avoids ambiguous
dimension semantics and allows the normal AgrocosmData reader to use an
explicit single-entry band mapping.

## 4. Validate before transfer

Check headers and sizes on the server:

```bash
ncdump -h /output/path/landuse_wheat_rainfed.nc
ncdump -h /output/path/temp_first_two_years.nc
```

Expected properties:

- management CFT/band dimension: `1`;
- longitude and latitude dimensions: unchanged;
- management time dimension: unchanged;
- climate time dimension: exactly 730 rows;
- all four climate files: identical time coordinate;
- CO₂ years: identical to the selected climate years.

Then point a local AgrocosmData catalog at the extracted files. For each
single-band management dataset use:

```toml
rainfed_bands = [1]
irrigated_bands = [1]
```

The second mapping is only a schema placeholder for these rainfed-only test
files; the production test must call readers with `irrigated = false`.

## 5. Scope

This extraction reduces time and CFT volume, not space. It is deliberately a
global fixture for testing the real `720 × 280` alignment, land-use mask,
compact cell ordering, streamed forcing, CPU/GPU execution, and output
reconstruction. It is not a scientifically complete historical experiment.

## 6. CPU production workflow

Copy and edit the production configuration:

```bash
cp examples/scripts/global_wheat_cpu.example.toml global_wheat_cpu.toml

julia --project=. examples/scripts/run_global_wheat_cpu.jl \
  /absolute/path/global_wheat_cpu.toml
```

The default `management.mode = "fixed"` repeats `fixed_year = 2015` throughout
the run. With `management.mode = "transient"`, each simulation year reads its
corresponding management row. Years before the file begins repeat its first
row, and years after it ends repeat its last row, matching LPJmL's boundary
behavior. Annual PHU is installed only when a new crop is sown, so a winter
crop already growing across 1 January retains the PHU assigned in its sowing
year. `landfrac > 0` selects crop cells but never scales single-cell processes.

For the global production experiment, point the `[climate]` file entries at
the complete 1901--2019 forcing files. Warm-up repeatedly uses
`warmup_climate_start_year = 1901` through `warmup_climate_end_year = 1930`,
while production still reads 2015--2016. Convergence compares soil state at
the same position in the 30-year forcing cycle (`t` versus `t-30`); it is not
evaluated from unlike adjacent climate years. Management remains fixed at the
2015 level during warm-up. The production configuration runs five complete
cycles, for a fixed total of 150 warm-up years.

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
  examples/scripts/run_global_wheat_gpu.jl \
  /absolute/path/global_wheat_gpu.toml
```

Set `run.device_id = 0` unless the scheduler exposes a different CUDA device.
Review `recommended_device_peak_gib` in `memory_preflight.toml` before the full
run. Validate 10 cells and then a bounded 1000-cell run before setting
`cell_limit = 0`.
