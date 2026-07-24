# Crop plant-available-water stress: the soil-moisture limiting factor β ∈ [0,1] used by crop
# photosynthesis and stomatal conductance. β is the fraction of plant-available water remaining
# between the wilting point and field capacity. The depth-integrated, root-weighted coupling to the
# Terrarium soil hydraulics is wired in the crop vegetation model; this is the tested
# scalar physics.

"""
    $(TYPEDSIGNATURES)

Soil-moisture limiting factor `β ∈ [0,1]` from the volumetric water content `θ`, the wilting-point
water content `θ_wilting`, and the field-capacity water content `θ_field_capacity`:
`β = clamp((θ − θ_wilting) / (θ_field_capacity − θ_wilting), 0, 1)`.
"""
@inline function soil_moisture_limiting_factor(θ::NF, θ_wilting::NF, θ_field_capacity::NF) where {NF}
    available = θ_field_capacity - θ_wilting
    β = (θ - θ_wilting) / max(available, eps(NF))
    return clamp(β, zero(NF), one(NF))
end

"""
    $(TYPEDSIGNATURES)

Plant-available water (m) in a soil layer of thickness `Δz` (m): the water between the wilting point
and the current content, floored at zero.
"""
@inline function plant_available_water(θ::NF, θ_wilting::NF, Δz::NF) where {NF}
    return max(zero(NF), (θ - θ_wilting)) * Δz
end
