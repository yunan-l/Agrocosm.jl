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
pft_ids = "all"

# One or both water systems.
water_systems = ["rainfed", "irrigated"]
```

```bash
julia --project=. examples/scripts/run_global_cfts_cpu.jl \
  /absolute/path/global_crops.toml
```

For CUDA, use `run_global_cfts_gpu.jl` with the same configuration.

The extracted management files are ordered as rainfed CFTs `1:12`, irrigated
CFTs `13:24`, and 12 residue bands shared by the two water systems. Before a
run, the selected land fractions are checked so that their sum cannot exceed
one in any geographic cell.

## Outputs and warm-up

Each batch writes below `batches/cft_XX_rainfed` or
`batches/cft_XX_irrigated`. The combined file
`global_cft_yield_START_END.nc` contains:

```text
yield(longitude, latitude, band, time)
pft_id(band)
irrigated(band)  # 0 = rainfed; 1 = irrigated
```

`yield` is deliberately unweighted: it is a per-patch yield and is not summed
or multiplied by `landfrac`. Full irrigation resets rooted soil layers to field
capacity daily, so water-balance closure is reported only for rainfed batches.

Each batch writes `warmup_soil_pool_allocation.nc`. Its allocation fields have
dimensions `(layer, cell, patch)` and its singleton `patch` coordinate carries
`pft_id` and `irrigated`; these metadata identify the CFT and water system that
calibrated the file. A soil-pool allocation calibrated for one patch must not
be reused for another CFT or water system.
