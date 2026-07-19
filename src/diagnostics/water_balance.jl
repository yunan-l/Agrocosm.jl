"""
    WaterBalance

Optional daily water-budget diagnostics. All fields have shape
`(number_of_days, number_of_cells)` and use millimetres of water.

`residual` is positive when water enters the model but is neither retained in
soil, snow, or surface-litter storage nor represented by a recorded outgoing
flux. Surface-litter interception is an internal transfer; litter evaporation
is an outgoing flux.
"""
mutable struct WaterBalance{M <: AbstractArray{<:AbstractFloat}}
    precipitation::M
    rain_after_snow::M
    soil_storage_before::M
    soil_storage_after::M
    soil_ice_storage_before::M
    soil_ice_storage_after::M
    snow_storage_before::M
    snow_storage_after::M
    litter_storage_before::M
    litter_storage_after::M
    snowmelt::M
    snow_sublimation::M
    snow_runoff::M
    unaccounted_snow_flux::M
    interception::M
    litter_interception::M
    litter_evaporation::M
    transpiration::M
    evaporation::M
    surface_runoff::M
    lateral_runoff::M
    bottom_drainage::M
    remaining_infiltration::M
    residual::M
end

"""
    init_water_balance(number_of_days, number_of_cells, device=identity; T=Float32)

Allocate an optional daily water-balance ledger on the selected device.
Currently the ledger supports rainfed simulations only.
"""
function init_water_balance(number_of_days::Integer,
                            number_of_cells::Integer,
                            device = identity;
                            T::Type{<:AbstractFloat} = Float32)
    allocate() = device(zeros(T, number_of_days, number_of_cells))

    return WaterBalance(
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(),
    )
end

function record_water_balance_start!(water_balance::WaterBalance,
                                     day_index::Integer,
                                     soil::Soil,
                                     precipitation)
    @views water_balance.precipitation[day_index, :] .= precipitation
    @views water_balance.soil_storage_before[day_index, :] .= vec(sum(soil.water.storage; dims = 1))
    @views water_balance.soil_ice_storage_before[day_index, :] .=
        vec(sum(soil.water.ice_storage; dims = 1))
    @views water_balance.snow_storage_before[day_index, :] .= soil.snow.pack
    @views water_balance.litter_storage_before[day_index, :] .=
        soil.surface_litter.water_storage
    return nothing
end

function record_water_balance_after_snow!(water_balance::WaterBalance,
                                          day_index::Integer,
                                          precipitation)
    @views water_balance.rain_after_snow[day_index, :] .= precipitation
    return nothing
end

function record_water_balance_end!(water_balance::WaterBalance,
                                   day_index::Integer,
                                   soil::Soil,
                                   crop::Crop)
    @views begin
        water_balance.soil_storage_after[day_index, :] .= vec(sum(soil.water.storage; dims = 1))
        water_balance.soil_ice_storage_after[day_index, :] .=
            vec(sum(soil.water.ice_storage; dims = 1))
        water_balance.snow_storage_after[day_index, :] .= soil.snow.pack
        water_balance.litter_storage_after[day_index, :] .=
            soil.surface_litter.water_storage
        water_balance.snowmelt[day_index, :] .= soil.snow.melt
        water_balance.snow_sublimation[day_index, :] .= soil.snow.sublimation
        water_balance.snow_runoff[day_index, :] .= soil.snow.runoff
        water_balance.interception[day_index, :] .= crop.water.interception
        water_balance.litter_interception[day_index, :] .=
            soil.surface_litter.interception
        water_balance.litter_evaporation[day_index, :] .=
            soil.surface_litter.evaporation
        water_balance.transpiration[day_index, :] .= vec(sum(crop.water.transpiration_layer; dims = 1))
        water_balance.evaporation[day_index, :] .= vec(sum(soil.water.evaporation; dims = 1))
        water_balance.surface_runoff[day_index, :] .= soil.water.surface_runoff
        water_balance.lateral_runoff[day_index, :] .= vec(sum(soil.water.lateral_runoff; dims = 1))
        water_balance.bottom_drainage[day_index, :] .= soil.water.bottom_drainage
        water_balance.remaining_infiltration[day_index, :] .= soil.water.infiltration

        water_balance.unaccounted_snow_flux[day_index, :] .=
            water_balance.snow_storage_before[day_index, :] .+
            water_balance.precipitation[day_index, :] .-
            water_balance.rain_after_snow[day_index, :] .-
            water_balance.snow_storage_after[day_index, :] .-
            water_balance.snow_sublimation[day_index, :] .-
            water_balance.snow_runoff[day_index, :]

        water_balance.residual[day_index, :] .=
            water_balance.soil_storage_before[day_index, :] .+
            water_balance.soil_ice_storage_before[day_index, :] .+
            water_balance.snow_storage_before[day_index, :] .+
            water_balance.litter_storage_before[day_index, :] .+
            water_balance.precipitation[day_index, :] .-
            water_balance.soil_storage_after[day_index, :] .-
            water_balance.soil_ice_storage_after[day_index, :] .-
            water_balance.snow_storage_after[day_index, :] .-
            water_balance.litter_storage_after[day_index, :] .-
            water_balance.interception[day_index, :] .-
            water_balance.litter_evaporation[day_index, :] .-
            water_balance.transpiration[day_index, :] .-
            water_balance.evaporation[day_index, :] .-
            water_balance.snow_sublimation[day_index, :] .-
            water_balance.snow_runoff[day_index, :] .-
            water_balance.surface_runoff[day_index, :] .-
            water_balance.lateral_runoff[day_index, :] .-
            water_balance.bottom_drainage[day_index, :]
    end

    return nothing
end
