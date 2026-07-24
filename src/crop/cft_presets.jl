# CFT presets: construct the crop processes from the LPJmL 12-CFT crop-trait registry (pft.jl). Each
# `PftParameters` entry (`crop_pft(id)` / `crop_pft(name)`) maps directly onto the crop process
# parameters — the phenological LAI trajectory, the C3/C4 pathway and temperature thresholds, and the
# base temperature — so a per-crop `CropVegetation` is `CropVegetation(NF, crop_pft(:maize))`.

"""$(TYPEDSIGNATURES) Crop phenology (LAI trajectory) for a CFT from its registry entry."""
function CropPhenology(::Type{NF}, pft::PftParameters) where {NF}
    return CropPhenology(
        NF;
        fphuc = NF(pft.fphuc),
        flaimaxc = NF(pft.flaimaxc),
        fphuk = NF(pft.fphuk),
        flaimaxk = NF(pft.flaimaxk),
        fphusen = NF(pft.fphusen),
        flaimaxharvest = NF(pft.flaimaxharvest),
        shapesenescencenorm = NF(pft.shapesenescencenorm),
        laimax = NF(pft.laimax),
    )
end

"""$(TYPEDSIGNATURES) Crop C3/C4 photosynthesis for a CFT (pathway + temperature thresholds)."""
function CropPhotosynthesis(::Type{NF}, pft::PftParameters) where {NF}
    pathway = pft.path == 1 ? C3Pathway() : C4Pathway()
    return CropPhotosynthesis(
        NF;
        pathway,
        T_CO2_low = NF(pft.temp_co2.low),
        T_CO2_high = NF(pft.temp_co2.high),
        T_photos_low = NF(pft.temp_photos.low),
        T_photos_high = NF(pft.temp_photos.high),
    )
end

"""$(TYPEDSIGNATURES) Crop heat-unit accumulation for a CFT (base temperature from the registry)."""
function CropPhenologyDynamics(::Type{NF}, pft::PftParameters) where {NF}
    # heat_unit_requirement is climate-derived (sowing→maturity heat sum) and kept at its default;
    # the base temperature comes from the CFT trait.
    return CropPhenologyDynamics(NF; base_temperature = NF(pft.basetemp.low))
end

"""
    $(TYPEDSIGNATURES)

Crop vegetation model for a CFT: assembles the per-crop phenology dynamics, LAI trajectory, and C3/C4
photosynthesis from the registry entry `pft` (e.g. `crop_pft(:maize)`).
"""
function CropVegetation(::Type{NF}, pft::PftParameters) where {NF}
    return CropVegetation(
        NF;
        phenology_dynamics = CropPhenologyDynamics(NF, pft),
        phenology = CropPhenology(NF, pft),
        photosynthesis = CropPhotosynthesis(NF, pft),
        # PFT specific leaf area is m²/gC; the carbon pool uses m²/kgC.
        carbon = CropCarbon(NF; specific_leaf_area = NF(pft.sla) * NF(1000)),
    )
end
