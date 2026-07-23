# Crop maintenance respiration (LPJmL crop `respiration`). Each living organ (root, storage, mobile
# pool) respires in proportion to its carbon, a respiration coefficient, its nitrogen:carbon ratio,
# and a Lloyd-Taylor soil/air-temperature response normalized to 1 at 10 °C. Total autotrophic
# respiration is this maintenance term plus growth respiration (see `growth_respiration.jl`). Tested
# scalar physics; assembled in the crop carbon coupling (plan Phase 5).

"""
    $(TYPEDEF)

Crop maintenance-respiration parameters.

Properties:
$(TYPEDFIELDS)
"""
@parameterized @kwdef struct CropMaintenanceRespiration{NF}
    "Organ maintenance-respiration coefficient"
    @param respcoeff::NF = 0.8 (bounds = Positive,)
    "Respiration scaling constant"
    @param k::NF = 0.0548 (bounds = Positive,)
    "Lloyd-Taylor activation parameter"
    @param e0::NF = 308.56 (bounds = Positive,)
    "Lloyd-Taylor reference-temperature offset"
    @param temp_response::NF = 56.02 (bounds = Positive,)
    "Root nitrogen:carbon ratio"
    @param nc_ratio_root::NF = 1 / 30 (bounds = Positive,)
    "Storage-organ nitrogen:carbon ratio"
    @param nc_ratio_storage::NF = 1 / 100 (bounds = Positive,)
    "Mobile-pool nitrogen:carbon ratio"
    @param nc_ratio_pool::NF = 1 / 100 (bounds = Positive,)
end

CropMaintenanceRespiration(::Type{NF}; kwargs...) where {NF} = CropMaintenanceRespiration{NF}(; kwargs...)

"""
    $(TYPEDSIGNATURES)

Lloyd-Taylor temperature response for maintenance respiration, normalized to 1 at 10 °C and gated to
zero at or below −15 °C.
"""
@inline function maintenance_temperature_response(m::CropMaintenanceRespiration{NF}, temperature::NF) where {NF}
    response = exp(m.e0 * (one(NF) / (m.temp_response + NF(10)) - one(NF) / (temperature + m.temp_response)))
    return ifelse(temperature ≥ NF(-15), response, zero(NF))
end

"""$(TYPEDSIGNATURES) Maintenance respiration of one organ: `carbon·respcoeff·k·nc_ratio·temperature_response`."""
@inline function organ_maintenance_respiration(m::CropMaintenanceRespiration{NF}, carbon::NF, nc_ratio::NF, temperature_response::NF) where {NF}
    return carbon * m.respcoeff * m.k * nc_ratio * temperature_response
end

"""
    $(TYPEDSIGNATURES)

Total maintenance respiration summed over the root (soil temperature), storage, and mobile-pool (air
temperature) organs, given their carbon pools and the soil/air temperatures (°C).
"""
@inline function maintenance_respiration(m::CropMaintenanceRespiration{NF}, root_carbon::NF, storage_carbon::NF, pool_carbon::NF, soil_temperature::NF, air_temperature::NF) where {NF}
    g_soil = maintenance_temperature_response(m, soil_temperature)
    g_air = maintenance_temperature_response(m, air_temperature)
    return organ_maintenance_respiration(m, root_carbon, m.nc_ratio_root, g_soil) +
        organ_maintenance_respiration(m, storage_carbon, m.nc_ratio_storage, g_air) +
        organ_maintenance_respiration(m, pool_carbon, m.nc_ratio_pool, g_air)
end
