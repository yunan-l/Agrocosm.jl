"""Daily diagnostics for soil enthalpy, ice, and freeze-thaw state."""
mutable struct ThermalBalance{M <: AbstractArray{<:AbstractFloat}}
    surface_energy_flux::M
    energy_residual::M
    column_energy::M
    total_ice_storage::M
    wilting_ice_storage::M
    available_ice_storage::M
    free_ice_storage::M
    ice_pool_residual::M
    maximum_frozen_fraction::M
    minimum_temperature::M
    maximum_temperature::M
end

function init_thermal_balance(number_of_days::Integer,
                              number_of_cells::Integer,
                              device = identity;
                              T::Type{<:AbstractFloat} = Float32)
    allocate() = device(zeros(T, number_of_days, number_of_cells))
    return ThermalBalance(
        allocate(), allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
    )
end

function record_thermal_balance!(thermal_balance::ThermalBalance,
                                 day_index::Integer,
                                 soil::Soil)
    @views begin
        thermal_balance.surface_energy_flux[day_index, :] .=
            soil.thermal.surface_energy_flux
        thermal_balance.energy_residual[day_index, :] .=
            soil.thermal.energy_residual
        thermal_balance.column_energy[day_index, :] .= vec(sum(
            soil.thermal.enthalpy .* reshape(
                soil.properties.layer_depth .* eltype(soil.properties.layer_depth)(0.001),
                :, 1,
            ); dims = 1,
        ))
        thermal_balance.total_ice_storage[day_index, :] .=
            vec(sum(soil.water.ice_storage; dims = 1))
        thermal_balance.wilting_ice_storage[day_index, :] .= vec(sum(
            soil.water.wilting_ice_fraction .* soil.water.wilting_storage;
            dims = 1,
        ))
        thermal_balance.available_ice_storage[day_index, :] .=
            vec(sum(soil.water.available_ice_storage; dims = 1))
        thermal_balance.free_ice_storage[day_index, :] .=
            vec(sum(soil.water.free_ice_storage; dims = 1))
        component_ice =
            soil.water.wilting_ice_fraction .* soil.water.wilting_storage .+
            soil.water.available_ice_storage .+
            soil.water.free_ice_storage
        thermal_balance.ice_pool_residual[day_index, :] .=
            vec(maximum(abs.(soil.water.ice_storage .- component_ice); dims = 1))
        thermal_balance.maximum_frozen_fraction[day_index, :] .=
            vec(maximum(soil.thermal.frozen_fraction; dims = 1))
        thermal_balance.minimum_temperature[day_index, :] .=
            vec(minimum(soil.thermal.temperature; dims = 1))
        thermal_balance.maximum_temperature[day_index, :] .=
            vec(maximum(soil.thermal.temperature; dims = 1))
    end
    return nothing
end
