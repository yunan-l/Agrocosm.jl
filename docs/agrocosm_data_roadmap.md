# AgrocosmData.jl roadmap

`AgrocosmData.jl` will provide the reproducible data layer for site, regional,
and global Agrocosm simulations. It will convert heterogeneous source files
into backend-neutral, cell-indexed arrays while keeping scientific processes,
spin-up dynamics, and device execution in `Agrocosm.jl`.

## Design decisions

- `grid.nc` is the only canonical spatial reference. Every record is keyed by
  its `cellid`; array order and source dimension names are never assumed.
- The current global rainfed-wheat experiment selects cells only where 2015
  `landfrac > 0`. The fraction is selection/provenance data, not a multiplier
  in model processes. All management values are fixed at 2015 levels.
- All allocated crop cells run together on one CPU or GPU when memory permits.
  Climate forcing is streamed in time blocks. Spatial batches are an explicit
  out-of-memory fallback, not the default execution model.
- The model and all management readers share LPJmL's 12-crop CFT order. Each
  dataset declares its own rainfed and irrigated band positions; total band
  count is never treated as the crop registry. `countries` is not an input.
- HWSD 2.x supplies layer-resolved SOC and total-N targets. Bulk density and
  coarse-fragment data are used only to convert concentrations to stocks.
  Existing `soilcode` lookup tables continue to supply Phase-2 hydraulic and
  thermal properties.
- The data package performs no GPU computation and no ecosystem spin-up. It
  returns CPU arrays plus metadata; `Agrocosm.jl` owns device transfer,
  hydraulic initialization, C/N pool partitioning, and state evolution.

## Current status

The package core is substantially complete. Milestones 1–5 have implemented
and fixture-tested contracts. A configuration-driven utility now extracts the
rainfed-wheat band for 2015 from the 64/32/24/16-band server files and the
2015–2016 daily forcing without loading the full datasets. Remaining work is
concentrated in full-grid HWSD quality control and the global smoke test. Large production files
remain on the server; their mappings are explicit code contracts tested with
small dimension-permuted fixtures.

## Milestone 1 — package and data contracts

Status: complete.

Create `AgrocosmData.jl` as a sibling Julia package with its own tests and a
small checked-in fixture dataset.

Deliverables:

- `DatasetCatalog` for configured server/local paths without hard-coded paths;
- versioned schemas for grid, static soil, management, climate, and restart
  data;
- backend-neutral batch objects compatible with the existing
  `InitialDataLoader` and `ClimateDataLoader` boundaries;
- source version, units, calendar, missing-value policy, and preprocessing
  provenance in every generated dataset.

Acceptance:

- the package loads without CUDA;
- the ten-cell example can be reconstructed through the new contracts without
  changing model results;
- dimension permutations and invalid/missing metadata fail with clear errors.

## Milestone 2 — canonical grid and crop masks

Status: complete in AgrocosmData. For the current experiment, the 2015 positive
fraction mask defines the fixed compact selection. Fraction and activity arrays
remain data products but are not model process inputs.

Build a reusable compact index from the `720 × 280` grid:

```text
compact_index ↔ cellid ↔ (longitude_index, latitude_index) ↔ (lon, lat)
```

Deliverables:

- exact coordinate validation for already aligned 0.5-degree inputs;
- conversion between gridded fields and compact `cell` arrays;
- `allocation_mask(landuse, cft, years)` and crop-fraction metadata;
- the 12-crop CFT registry plus explicit rainfed/irrigated band maps for every
  management file.

Acceptance:

- compact-to-grid round trips preserve every value and `cellid`;
- masks exclude missing grid cells and include every cell with positive
  2015 land-use fraction for the current production experiment;
- changing NetCDF dimension names or order does not change the result.

## Milestone 3 — current soil and management inputs

Status: complete at the reader, mapping, and fixture-equivalence level. The
production files do not need to be copied into the repository or materialized
for every CFT locally.

Move the existing NetCDF reading and lookup work out of the model package.

Deliverables:

- readers for `soilcode`, soil pH, land use, sowing state, PHU, fertilizer,
  manure, and residue fraction;
- model configuration for globally enabling or disabling tillage, with no
  `with_tillage` dataset dependency;
- fertilizer input required only for `fertilizer = :yes`; `:no` and `:auto`
  do not read it, while manure remains independently configured;
- the existing soil-code-to-property lookup represented as a versioned
  Agrocosm dataset rather than executable LPJmL input;
- unit and range checks before arrays reach model initialization;
- crop-specific reads that select one CFT without materializing every CFT.

Acceptance:

- the ten-cell baseline exactly matches the current loader;
- all static inputs align to the canonical `cellid` sequence;
- every climate, management, soil, and HWSD lookup obtains coordinates from
  the same canonical `grid.nc`; no process reconstructs the grid origin;
- no `countries` file is opened or required.

## Milestone 4 — HWSD C/N preprocessing

Implement a one-time offline HWSD 2.x conversion pipeline.

Status: the versioned numerical core is implemented. It converts HWSD v2.x
organic-carbon, total-nitrogen, bulk-density, and coarse-fragment attributes
to layer stocks; maps the official seven layers to the five Agrocosm layers;
and incrementally aggregates raster tiles by spherical pixel area. The
unsupported 2–3 m interval extends the 1.5–2 m stock density by default and is
always marked uncertain. Compact NetCDF target writing, coverage maps, and
conservation totals are included. The official HWSD v2.01 mapping-unit
database/raster adapter is implemented and validated on real cells. Incomplete
cells use a bounded nearest-complete-profile fallback with donor provenance.
Component mixing now preserves original `SHARE`, treats layers below
`ROOT_DEPTH` as structural zero, and excludes NODATA, water, and glacier from
the soil-area denominator. Generating the complete canonical-grid product and
reviewing its global missing-data and conservation summaries remain production
tasks.

Deliverables:

- area-conservative aggregation from HWSD resolution to the Agrocosm
  0.5-degree grid;
- vertical overlap mapping from the seven HWSD layers to the Agrocosm layers
  `0–0.2`, `0.2–0.5`, `0.5–1`, `1–2`, and `2–3 m`;
- layer SOC and total-N stocks in `gC m⁻²` and `gN m⁻²`;
- an explicit, versioned rule and uncertainty flag for the unsupported
  `2–3 m` interval;
- global conservation summaries, missing-data maps, and point comparisons for
  validation.

The preprocessing output contains total SOC/N targets, not `fast`/`slow` pools
or litter. Agrocosm's `soil_initial_state` applies the explicit interim 40:60
fast/slow initialization and zero litter before agricultural warm-up.

Acceptance:

- horizontal and vertical aggregation closes against the source stocks within
  documented numerical tolerance;
- source versions and every conversion are recoverable from metadata;
- rerunning preprocessing produces byte-stable metadata and numerically
  identical arrays.

## Milestone 5 — global climate streaming

Replace eager whole-file reads with time-blocked access over the fixed active
cell set.

Status: complete for the currently required forcing contract: daily
temperature, precipitation, net longwave, and downward shortwave, with annual
global CO₂ matched to each block. Gregorian leap days are normalized to the
model's 365-day calendar, forcing units are converted to the canonical model
contract, and one-block threaded prefetch is implemented. Real server-I/O
benchmarking remains part of production hardening.

Deliverables:

- climate readers for temperature, precipitation, shortwave, longwave, and CO₂
  with calendar and unit normalization; add wind/humidity only when a selected
  process configuration declares them required;
- configurable monthly, annual, or multi-year time blocks;
- optional canonical `time × cell` caches where benchmarks show that direct
  source reads are too slow;
- reusable host buffers and a prefetch interface suitable for overlapping CPU
  I/O with GPU execution.

Acceptance:

- concatenated blocks are numerically identical to an eager read;
- state and climate-buffer continuity are unchanged across block boundaries;
- memory use is bounded by active model state plus one or two forcing blocks,
  not total simulation duration.

## Milestone 6 — global runner and restart workflow

Integrate the data package with the Agrocosm simulation API.

Default execution:

```text
grid + CFT + years
  → cells with 2015 landfrac > 0
  → all active cells on one backend
  → streamed climate blocks
  → cellid-keyed output/checkpoint
```

Deliverables:

- memory estimation and automatic selection between whole-mask execution and
  spatial fallback batches;
- fixed 2015 management reused across production years;
- output reconstruction to `720 × 280` using `cellid`;
- per-shard checkpoints and deterministic merge for fallback batching;
- reproducibility metadata containing CFT, years, masks, source versions, and
  model/data schema versions.

Acceptance:

- whole-mask and spatially batched executions are equal after reassembly;
- CPU and GPU runs use the same compact cell ordering;
- restart across a climate-block boundary matches an uninterrupted run;
- a one-year global crop smoke test completes with bounded memory and closed
  model balance diagnostics.

Status: production validation. `estimate_memory` reports model, diagnostic,
output, scratch, forcing-transfer, warm-up history, cached forcing, and
prefetched-host memory. Streamed selected output, identity-checked
checkpoint/restart, asynchronous one-block prefetch, canonical HWSD
preprocessing/QC, and full selected-domain CPU/GPU runners are implemented.
The first ten selected cells pass the HWSD-backed CPU workflow. Non-interactive
Slurm templates now submit full CPU and single-GPU jobs. Current remaining work
is server backend validation and, only if whole-mask device memory is
insufficient, deterministic spatial fallback batching.

## Milestone 7 — spin-up handoff and production hardening

Connect HWSD targets and streamed forcing to the Agrocosm spin-up workflow
without moving scientific state evolution into the data package.

Interim direct-initialization status: equilibrium spin-up is deferred.
`soil_initial_state`
conservatively partitions HWSD SOC and total N into the current fast/slow
pools, initializes litter to zero, accounts for the default initial mineral-N
pools, and initializes water at field capacity. New data use the neutral
top-level `initial_state` contract.
compatibility path only.

`agricultural_warmup!` now cycles the configured 1901–1930 forcing for five
complete cycles (150 years) without advancing the production
clock or retaining production outputs. It supports target-constrained annual
corrections, same-climate-phase convergence, a strict production checkpoint
gate, restartable climate-block readers, and annual C/N/water diagnostics. The
150-year global result will determine whether a reusable target-constrained
pool-allocation product is needed.

Deliverables:

- loading and writing native Agrocosm spin-up checkpoints;
- HWSD SOC/N targets, field-capacity water initialization, and convergence
  diagnostics exposed through a stable handoff contract;
- no runtime dependency on a legacy external initial-state schema
  and restart tests pass;
- server deployment documentation, data catalog examples, and performance
  benchmarks for representative CFT masks.

Acceptance:

- `swc`, `litc`, `fastc`, `slowc`, `litn`, `fastn`, and `slown` can all be
  supplied from a native Agrocosm checkpoint;
- expanding to new grid cells requires only source-data preprocessing and
  Agrocosm spin-up, never an LPJmL run;
- the legacy ten-cell input remains available solely as a regression fixture.

## Deferred work

- simultaneous CPU/GPU work-queue scheduling;
- multi-GPU and distributed-memory execution;
- alternative soil hydraulic datasets replacing the current soil-code lookup;
- dynamic compaction that changes allocated cells during a simulation.

These should follow the single-backend global runner because they do not alter
the canonical data or cell-index contracts.
