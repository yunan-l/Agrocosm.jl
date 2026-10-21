function _enzyme_evaporation!(
    state::Agrocosm.ModelState,
    pet_eeq,
    lpjmlparams,
    layer_depth,
)
    T = eltype(pet_eeq)
    crop_fpar = Agrocosm.crop_canopy_auxiliary(state).fpar
    crop_transpiration = Agrocosm.crop_fluxes(state).water.transpiration_layer
    crop_canopy_wet = Agrocosm.crop_canopy_auxiliary(state).canopy_wet
    soil_water = Agrocosm.soil_water_auxiliary(state)
    soil_water_state = Agrocosm.soil_water_prognostic(state)
    soil_water_flux = Agrocosm.soil_water_fluxes(state)
    surface_litter = Agrocosm.soil_surface_litter_prognostic(state)
    surface_litter_aux = Agrocosm.soil_surface_litter_auxiliary(state)
    surface_litter_flux = Agrocosm.soil_surface_litter_fluxes(state)
    priestley_taylor = lpjmlparams.PRIESTLEY_TAYLOR
    for cell in eachindex(pet_eeq)
        evap_energy = pet_eeq[cell] * priestley_taylor *
            max(one(T) - crop_fpar[cell], T(0.1))
        crop_transpiration_sum = crop_transpiration[1, cell] +
            crop_transpiration[2, cell] + crop_transpiration[3, cell] +
            crop_transpiration[4, cell] + crop_transpiration[5, cell]
        available_evaporation = max(
            pet_eeq[cell] * priestley_taylor *
                (one(T) - crop_canopy_wet[cell]) - crop_transpiration_sum,
            zero(T),
        )
        surface_litter_flux.evaporation[cell] = Agrocosm.compute_litter_evaporation(
            evap_energy,
            available_evaporation,
            surface_litter.water_storage[cell],
            surface_litter_aux.water_capacity[cell],
            surface_litter.cover[cell],
        )
        surface_litter.water_storage[cell] -= surface_litter_flux.evaporation[cell]

        liquid_1 = soil_water.relative_content[1, cell] *
            soil_water.holding_capacity_storage[1, cell] + soil_water.free_water[1, cell]
        liquid_2 = soil_water.relative_content[2, cell] *
            soil_water.holding_capacity_storage[2, cell] + soil_water.free_water[2, cell]
        liquid_3 = soil_water.relative_content[3, cell] *
            soil_water.holding_capacity_storage[3, cell] + soil_water.free_water[3, cell]
        liquid_4 = soil_water.relative_content[4, cell] *
            soil_water.holding_capacity_storage[4, cell] + soil_water.free_water[4, cell]
        liquid_5 = soil_water.relative_content[5, cell] *
            soil_water.holding_capacity_storage[5, cell] + soil_water.free_water[5, cell]
        evaporating_liquid_1 = max(liquid_1 - crop_transpiration[1, cell], zero(T))
        evaporating_liquid_2 = max(liquid_2 - crop_transpiration[2, cell], zero(T))
        evaporating_liquid_3 = max(liquid_3 - crop_transpiration[3, cell], zero(T))
        evaporating_liquid_4 = max(liquid_4 - crop_transpiration[4, cell], zero(T))
        evaporating_liquid_5 = max(liquid_5 - crop_transpiration[5, cell], zero(T))
        evap_ratio = zero(T)
        if evap_energy > T(1e-5) && available_evaporation > T(1e-5)
            soil_depth_evap = T(lpjmlparams.soildepth_evap)
            fraction_1 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[1])
            tmpwater = evaporating_liquid_1 * fraction_1
            water_evap = (
                evaporating_liquid_1 + soil_water_state.available_ice_storage[1, cell] +
                soil_water_state.free_ice_storage[1, cell]
            ) * fraction_1
            holding_evap = soil_water.holding_capacity_storage[1, cell] * fraction_1
            soil_depth_evap -= layer_depth[1]
            fraction_2 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[2])
            tmpwater += evaporating_liquid_2 * fraction_2
            water_evap += (
                evaporating_liquid_2 + soil_water_state.available_ice_storage[2, cell] +
                soil_water_state.free_ice_storage[2, cell]
            ) * fraction_2
            holding_evap += soil_water.holding_capacity_storage[2, cell] * fraction_2
            soil_depth_evap -= layer_depth[2]
            fraction_3 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[3])
            tmpwater += evaporating_liquid_3 * fraction_3
            water_evap += (
                evaporating_liquid_3 + soil_water_state.available_ice_storage[3, cell] +
                soil_water_state.free_ice_storage[3, cell]
            ) * fraction_3
            holding_evap += soil_water.holding_capacity_storage[3, cell] * fraction_3
            soil_depth_evap -= layer_depth[3]
            fraction_4 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[4])
            tmpwater += evaporating_liquid_4 * fraction_4
            water_evap += (
                evaporating_liquid_4 + soil_water_state.available_ice_storage[4, cell] +
                soil_water_state.free_ice_storage[4, cell]
            ) * fraction_4
            holding_evap += soil_water.holding_capacity_storage[4, cell] * fraction_4
            soil_depth_evap -= layer_depth[4]
            fraction_5 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[5])
            tmpwater += evaporating_liquid_5 * fraction_5
            water_evap += (
                evaporating_liquid_5 + soil_water_state.available_ice_storage[5, cell] +
                soil_water_state.free_ice_storage[5, cell]
            ) * fraction_5
            holding_evap += soil_water.holding_capacity_storage[5, cell] * fraction_5
            evap_ratio = Agrocosm.compute_lpjml_managed_soil_evaporation_ratio(
                evap_energy,
                available_evaporation,
                water_evap,
                tmpwater,
                holding_evap,
                surface_litter.cover[cell],
            )
        end

        soil_depth_evap = T(lpjmlparams.soildepth_evap)
        fraction_1 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[1])
        soil_water_flux.evaporation[1, cell] = evaporating_liquid_1 * evap_ratio * fraction_1
        soil_depth_evap -= layer_depth[1]
        fraction_2 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[2])
        soil_water_flux.evaporation[2, cell] = evaporating_liquid_2 * evap_ratio * fraction_2
        soil_depth_evap -= layer_depth[2]
        fraction_3 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[3])
        soil_water_flux.evaporation[3, cell] = evaporating_liquid_3 * evap_ratio * fraction_3
        soil_depth_evap -= layer_depth[3]
        fraction_4 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[4])
        soil_water_flux.evaporation[4, cell] = evaporating_liquid_4 * evap_ratio * fraction_4
        soil_depth_evap -= layer_depth[4]
        fraction_5 = min(one(T), max(soil_depth_evap, zero(T)) / layer_depth[5])
        soil_water_flux.evaporation[5, cell] = evaporating_liquid_5 * evap_ratio * fraction_5
    end
    return nothing
end

function _enzyme_partition_soil_water_ice!(state::Agrocosm.ModelState)
    water = Agrocosm.soil_water_prognostic(state)
    auxiliary = Agrocosm.soil_water_auxiliary(state)
    for cell in eachindex(Agrocosm.crop_fluxes(state).carbon.gross_assimilation)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 1, cell)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 2, cell)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 3, cell)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 4, cell)
        _enzyme_partition_soil_water_ice_layer!(water, auxiliary, 5, cell)
    end
    return nothing
end

function _enzyme_partition_soil_water_ice_layer!(water, auxiliary, layer, cell)
    T = eltype(water.storage)
    total = max(water.storage[layer, cell] + water.ice_storage[layer, cell], zero(T))
    ice = clamp(water.ice_storage[layer, cell], zero(T), total)
    pwp_fraction, available_ice, free_ice, relative_water, free_water =
        Agrocosm.lpjml_water_ice_partition(
            total,
            ice,
            auxiliary.wilting_storage[layer, cell],
            auxiliary.holding_capacity_storage[layer, cell],
        )
    water.ice_storage[layer, cell] = ice
    water.storage[layer, cell] = total - ice
    water.wilting_ice_fraction[layer, cell] = pwp_fraction
    water.available_ice_storage[layer, cell] = available_ice
    water.free_ice_storage[layer, cell] = free_ice
    auxiliary.relative_content[layer, cell] = relative_water
    auxiliary.free_water[layer, cell] = free_water
    return nothing
end

function _enzyme_soil_evapotranspiration!(
    state::Agrocosm.ModelState;
    irrigation::Bool = false,
)
    water = Agrocosm.soil_water_prognostic(state)
    storage = water.storage
    transpiration = Agrocosm.crop_fluxes(state).water.transpiration_layer
    evaporation = Agrocosm.soil_water_fluxes(state).evaporation
    if irrigation
        field_capacity = Agrocosm.soil_water_auxiliary(state).field_capacity
        layer_depth = Agrocosm.soil_properties(state).layer_depth
        ice_storage = water.ice_storage
        for cell in axes(storage, 2), layer in axes(storage, 1)
            target_storage = field_capacity[layer, cell] * layer_depth[layer]
            storage[layer, cell] = max(
                target_storage - ice_storage[layer, cell],
                zero(eltype(storage)),
            )
        end
    else
        for cell in eachindex(Agrocosm.crop_fluxes(state).carbon.gross_assimilation)
            storage[1, cell] -= transpiration[1, cell] + evaporation[1, cell]
            storage[2, cell] -= transpiration[2, cell] + evaporation[2, cell]
            storage[3, cell] -= transpiration[3, cell] + evaporation[3, cell]
            storage[4, cell] -= transpiration[4, cell] + evaporation[4, cell]
            storage[5, cell] -= transpiration[5, cell] + evaporation[5, cell]
        end
    end
    _enzyme_partition_soil_water_ice!(state)
    return nothing
end

"""Keep the AD root profile linked to the active CFT beta_root parameter.

The production initializer evaluates this profile before the daily transition,
which is correct for ordinary runs but makes `beta_root` invisible when the
post-sowing state is treated as fixed. Reassigning the same profile inside the
opt-in AD path preserves baseline primal parity and exposes the intended
initial-condition sensitivity without changing production initialization.
"""
function _enzyme_apply_root_distribution!(state::Agrocosm.ModelState, beta_root)
    distribution = state.inputs.crop.root.distribution
    T = eltype(distribution)
    beta = T(beta_root)
    beta20 = beta ^ 20
    beta50 = beta ^ 50
    beta100 = beta ^ 100
    beta200 = beta ^ 200
    beta300 = beta ^ 300
    total = one(T) - beta300
    distribution[1] = (one(T) - beta20) / total
    distribution[2] = (beta20 - beta50) / total
    distribution[3] = (beta50 - beta100) / total
    distribution[4] = (beta100 - beta200) / total
    distribution[5] = (beta200 - beta300) / total
    return nothing
end

"""Keep same-precision AD parameter bundles on the differentiable path."""
@inline function _enzyme_model_parameters(
    ::Type{T},
    parameters::Agrocosm.ModelParameters{T},
) where {T <: AbstractFloat}
    return parameters
end

@inline function _enzyme_model_parameters(
    ::Type{T},
    parameters::Agrocosm.ModelParameters,
) where {T <: AbstractFloat}
    return Agrocosm.convert_precision(T, parameters)
end
