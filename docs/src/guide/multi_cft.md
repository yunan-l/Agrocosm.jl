# Multi-CFT simulations

The multi-CFT drivers run any selected subset of the 12 LPJmL-compatible crop
functional types (CFTs) under rainfed, full irrigation, or both. Each selected
CFT and water system is an independent crop--soil patch. A geographic grid cell
may therefore occur in several output bands, each with its own crop, soil,
litter, management, warm-up, and checkpoint state.

The first implementation groups identical CFT/water-system patches into
separate batches and reuses the established single-crop process kernels.
Climate, coordinates, soil code, and pH are read once per batch; they do not
imply a shared soil state between crop patches.

## Configuration

The single-crop runners use the standard 24-band management products and, by
default, simulate rainfed CFT 1 from their first band. Add a `[cfts]` section
and use the multi-CFT driver to select several patches:

```toml
[cfts]
# A single CFT, an arbitrary unique subset, or "all".
cft_ids = "all"

# One or both water systems.
water_systems = ["rainfed", "irrigated"]
```

To calibrate and initialize a separate soil-pool state for every selected patch
in one workflow, add a `[free_warmup]` section. The supplied
`global_cfts_allocation_validation.example.toml` is the canonical template.
Every calibration phase runs a fixed 600-year target-constrained warm-up and
writes an allocation. The free phase reloads that allocation, runs for at least
600 years without target correction, and stops by 1500 years only after its
configured convergence fraction is met. Production then restores the final
free warm-up checkpoint and runs the requested years with checkpoint/restart.

```bash
julia --project=. examples/scripts/run_global_cfts_cpu.jl \
  /absolute/path/global_crops.toml
```

For CUDA, use `run_global_cfts_gpu.jl` with the same configuration.
Both drivers execute the selected homogeneous CFT/water-system batches
sequentially. This preserves backend-equivalent process execution during the
current validation stage; it is not yet a fused compact multi-patch GPU domain.

The extracted management files are ordered as rainfed CFTs `1:12`, irrigated
CFTs `13:24`, and 12 residue bands shared by the two water systems. Before a
run, the selected land fractions are checked so that their sum cannot exceed
one in any geographic cell.

## Outputs and warm-up

Each calibration/free-equilibrium batch writes separately below
`batches/cft_XX_rainfed/{calibration,production}` or
`batches/cft_XX_irrigated/{calibration,production}`. The combined file
`global_cft_yield_START_END.nc` contains:

```text
yield(longitude, latitude, band, time)
landfrac(longitude, latitude, band, time)
cft_id(band)
irrigated(band)  # 0 = rainfed; 1 = irrigated
landfrac_sum(longitude, latitude, time)
```

`yield` is deliberately unweighted: it is a per-patch yield and is not summed
or multiplied by `landfrac`; `landfrac_sum` is retained only as a coverage QC.
`cft_batch_manifest.toml` records the allocation and production paths for every
batch. Full irrigation resets rooted soil layers to field capacity daily, so
water-balance closure is reported only for rainfed batches.

Each batch writes `warmup_soil_pool_allocation.nc`. Its allocation fields have
dimensions `(layer, cell, patch)` and its singleton `patch` coordinate carries
`cft_id` and `irrigated`; these metadata identify the CFT and water system that
calibrated the file. A soil-pool allocation calibrated for one patch must not
be reused for another CFT or water system.
