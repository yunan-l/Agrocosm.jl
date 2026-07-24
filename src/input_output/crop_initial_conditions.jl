# Reading Agrocosm's legacy per-cell crop initial-condition / site file (the `initial_wheat.jld2`
# layout) and mapping the entries that carry cleanly onto the Terrarium crop model. The file stores,
# per grid cell, the crop management (sowing date, phenological-heat-unit requirement, residue
# fraction), the soil texture, and the LPJmL initial soil state (`u0`).

# Nominal LPJmL soil-layer thicknesses (m), used to convert the per-area carbon pools to densities.
const LPJML_LAYER_THICKNESS = (0.2, 0.3, 0.5, 1.0, 1.0)

"""
    $(TYPEDSIGNATURES)

Read Agrocosm's legacy initial-condition file (the `initial_wheat.jld2` layout) and return **per-cell
vectors** of the crop site setup mapped to the Terrarium crop model. Index the vectors by grid cell.

Fields:
- `sowing_day` — LPJmL sowing day of year (`crop.sdate`)
- `heat_unit_requirement` — phenological heat-unit requirement, °C·days (`crop.phu`), for
  [`CropPhenologyDynamics`](@ref) / `CropVegetation(NF, pft; heat_unit_requirement = …)`
- `residue_fraction` — fraction of residue returned to the soil at harvest (`crop.residuefrac`)
- `fast_carbon`, `slow_carbon` — initial fast/slow soil carbon density (kgC/m³) for
  [`CropSoilBiogeochemistry`](@ref), the column mean of the LPJmL per-layer pools converted with the
  nominal layer thicknesses
- `latitude` — cell latitude (degrees)

The file also stores soil texture and the LPJmL mineral-nitrogen pools; the mineral-N pools use LPJmL
units/scaling (and contain fill values) and are left to the caller.
"""
function load_crop_initial_conditions(path::AbstractString)
    data = JLD2.load(path, "initial_data")
    u0 = data.initialLPJmL.u0
    ncells = length(data.latitude)
    # Column-mean carbon density (kgC/m³) from the LPJmL per-layer per-area pools (gC/m²).
    function column_density(pool)
        nlayers = size(pool, 1)
        return [
            sum(Float64(pool[l, c]) / LPJML_LAYER_THICKNESS[l] / 1000 for l in 1:nlayers) / nlayers
                for c in 1:ncells
        ]
    end
    return (
        sowing_day = Int.(round.(Float64.(vec(data.crop.sdate[1, :])))),
        heat_unit_requirement = Float64.(data.crop.phu),
        residue_fraction = Float64.(data.crop.residuefrac),
        fast_carbon = column_density(u0.fastc),
        slow_carbon = column_density(u0.slowc),
        latitude = Float64.(data.latitude),
    )
end
