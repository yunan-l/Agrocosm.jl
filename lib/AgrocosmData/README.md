# AgrocosmData.jl

`AgrocosmData.jl` is the backend-neutral input-data layer for Agrocosm. It
aligns external NetCDF data to the canonical Agrocosm `cellid` grid, constructs
land-use/PFT selections, and loads compact soil and management arrays.

The package intentionally does not depend on CUDA or run scientific processes.
Agrocosm.jl remains responsible for device transfer, state initialization, and
simulation.

The implementation roadmap is maintained in
[`docs/agrocosm_data_roadmap.md`](../../docs/agrocosm_data_roadmap.md).

```julia
using AgrocosmData

catalog = load_catalog("catalog.toml")
grid = read_grid(dataset(catalog, :grid))

landuse = read_management(catalog, :landuse, grid, 1; years = 2000:2019)
crop_mask = build_crop_mask(grid, landuse.values)
soil = read_soil_data(catalog, grid; selection = crop_mask.selection)

reader = climate_blocks(
    catalog,
    grid;
    selection = crop_mask.selection,
    start_year = 2000,
    end_year = 2019,
    block_days = 31,
)
for block in reader
    forcing = climate_forcing(block) # bounded time × cell arrays
    # Transfer `forcing` to the selected backend and advance the model block.
end
```

`crop_mask.selection` has a stable, `cellid`-sorted order for the full run.
`crop_mask.active` retains the year-specific land-use state for switching crop
processes on and off without reallocating soil state.

The catalog uses the canonical 12-crop registry and declares separate band
positions for each file. Pass `irrigated=true` to select its irrigated band.
Mineral fertilizer follows the LPJmL-style `:no`, `:yes`, and `:auto` modes;
only `:yes` requires fertilizer input. Manure remains an independent switch and
requires manure input when enabled. Tillage is a model configuration switch,
not a data input.

Climate blocks read only the requested daily rows from `temp`, `prec`, `lwnet`,
and `swdown`. Annual global CO₂ from the two-column text file is matched by
calendar year and emitted as a small daily vector, so arbitrary block boundaries
remain correct.
