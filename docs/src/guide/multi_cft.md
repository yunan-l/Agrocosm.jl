# CFT simulations

The CFT drivers run one, several, or all 12 LPJmL-compatible crop functional
types (CFTs) under rainfed, full irrigation, or both. A configuration without a
`[cfts]` section defaults to rainfed CFT 1. Each selected CFT and water system
is an independent crop--soil patch. A geographic grid cell may therefore occur
in several output bands, each with its own crop, soil, litter, management,
warm-up, and checkpoint state.

The first implementation groups identical CFT/water-system patches into
separate batches and reuses the established single-crop process kernels.
Climate, coordinates, soil code, and pH are read once per batch; they do not
imply a shared soil state between crop patches.

## Configuration

Use `run_global_cfts_cpu.jl` or `run_global_cfts_gpu.jl` for both single- and
multi-CFT global testing. To request one crop explicitly, use `cft_ids = [1]`;
an omitted `[cfts]` section has the same rainfed-CFT-1 default. Add a `[cfts]`
section to select several patches:

```toml
[cfts]
# A single CFT, an arbitrary unique subset, or "all".
cft_ids = "all"

# One or both water systems.
water_systems = ["rainfed", "irrigated"]
```

`run_global_wheat_cpu.jl` remains the shared low-level batch executor used by
the CFT driver. It is retained for compatibility, but new server tests should
use the unified CFT entry points above.

To calibrate and initialize a separate soil-pool state for every selected patch
in one workflow, add a `[free_warmup]` section. The supplied
`global_cfts_allocation_validation.example.toml` is the canonical template.
Every calibration phase runs a fixed 600-year target-constrained warm-up and
writes an allocation. The default free-settling phase reloads that allocation,
runs five 30-year climate cycles (150 years) without target correction, then
uses the resulting state for production. It writes the convergence and drift
diagnostics but does not require free equilibrium; therefore its outputs are a
production baseline, not a claimed long-term soil equilibrium.

```bash
julia --project=. scripts/run_global_cfts_cpu.jl \
  /absolute/path/global_crops.toml
```

For CUDA, use `run_global_cfts_gpu.jl` with the same configuration.
Both drivers execute the selected homogeneous CFT/water-system batches
sequentially. This preserves backend-equivalent process execution during the
current validation stage; it is not yet a fused compact multi-patch GPU domain.

To run exactly one patch in a scheduler array or rerun one failed patch, pass
its CFT identifier and water system on the command line:

```bash
julia --project=. scripts/run_global_cfts_gpu.jl \
  /absolute/path/global_crops.toml 1 rainfed
```

With `resume_completed_batches = true` under `[run]`, a valid completed batch
is skipped. If calibration completed but the later free warm-up or production
stage stopped, the driver verifies the recorded allocation fingerprint and
resumes at production rather than recalibrating the 600-year allocation.

The extracted management files are ordered as rainfed CFTs `1:12`, irrigated
CFTs `13:24`, and 12 residue bands shared by the two water systems. Before a
run, the selected land fractions are checked so that their sum cannot exceed
one in any geographic cell.

## Outputs and warm-up

Each calibration/free-equilibrium batch writes separately below
`batches/cft_XX_rainfed/{calibration,production}` or
`batches/cft_XX_irrigated/{calibration,production}`. The production NetCDF in
each batch is the authoritative yield output; it is deliberately unweighted
and is not merged or multiplied by `landfrac`. `cft_batch_manifest.toml`
records the allocation, production directory, and production NetCDF path for
every batch. This keeps each crop--soil patch independent and avoids inventing
a cross-CFT aggregation rule at the model-output boundary.

`landfrac` is still read during input selection to identify active cells. It is
not used to rescale the process calculations or to construct a combined yield
file. Full irrigation resets rooted soil layers to field capacity daily, so
water-balance closure is reported only for rainfed batches.

Each batch writes `warmup_soil_pool_allocation.nc`. Its allocation fields have
dimensions `(layer, cell, patch)` and its singleton `patch` coordinate carries
`cft_id` and `irrigated`; these metadata identify the CFT and water system that
calibrated the file. A soil-pool allocation calibrated for one patch must not
be reused for another CFT or water system.
