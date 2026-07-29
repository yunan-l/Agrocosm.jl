# Multi-CFT simulations

The CPU multi-CFT runner simulates any selected subset of the 12
LPJmL-compatible crop functional types (CFTs) under rainfed, full irrigation,
or both water systems. Selecting all CFTs and both systems produces at most 24
independent crop-soil patches.

Each patch has independent crop, soil, litter, management, warm-up, and
checkpoint state. Climate, coordinates, soil-code, and pH inputs are shared by
the patches that occupy the same geographical grid cell. A grid cell can
therefore appear in more than one output band.

## Configuration

Start from `examples/scripts/global_wheat_cpu.example.toml` and configure the
requested crop systems:

```toml
[cfts]
# Any non-empty, unique subset of 1:12 is accepted.
pft_ids = [1, 2, 3]

# Select one or both modes.
water_systems = ["rainfed", "irrigated"]
```

Run the multi-CFT driver:

```bash
julia --project=. examples/scripts/run_global_cfts_cpu.jl \
  /absolute/path/global_crops.toml
```

The runner expects the extracted 2015 management files with the following
fixed ordering: bands `1:12` are rainfed CFTs, bands `13:24` are irrigated
CFTs, and the 12 residue bands are shared between the two water systems.

## Output contract

Each selected CFT and water system runs in a separate directory under
`batches/cft_XX_rainfed` or `batches/cft_XX_irrigated`. The final file is
`global_cft_yield_START_END.nc` with:

```text
yield(longitude, latitude, band, time)
pft_id(band)
irrigated(band)  # 0 = rainfed, 1 = irrigated
```

`yield` has `aggregation = none`: values are not weighted by `landfrac` and
are not summed across CFTs or grid cells. The selected patch land fractions
are checked before execution and must not exceed one per geographical cell.

## Full irrigation assumption

Full irrigation resets rooted soil layers to field capacity each day. It is an
unconstrained external water source, not a river, reservoir, withdrawal, or
delivery model. Accordingly, daily water-balance closure is evaluated for
rainfed patches only; irrigated batch reports explicitly record that this
diagnostic is not evaluated. Carbon, nitrogen, and thermal process paths are
otherwise identical.
