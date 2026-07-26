"""
harvest_crop!(crop, soil, output, residue_fraction, day)

Handle harvest-day biomass removal, residue transfer, and crop state reset.
"""
function harvest_crop!(crop,
                       soil,
                       output::Output,
                       residue_frac::AbstractArray{T},
                       day::Int;
                       output_row::Union{Nothing, Integer} = nothing,
                       annual_output_row::Union{Nothing, Integer} = nothing
) where {T <: AbstractFloat}
    annual_row = something(annual_output_row, 0)

    launch_1D!(
        harvest_state_kernel!,
        crop_events(crop).harvest,
        output.annual.harvest_date,
        crop_prognostic(crop).phenology.harvesting_previous,
        crop_prognostic(crop).phenology.harvesting,
        crop_prognostic(crop).phenology.is_growing,
        crop_fluxes(crop).carbon.yield,
        crop_fluxes(crop).carbon.harvest_export,
        output.annual.yield,
        crop_prognostic(crop).carbon.storage,
        crop_prognostic(crop).carbon.leaf,
        crop_prognostic(crop).carbon.pool,
        crop_prognostic(crop).carbon.root,
        crop_fluxes(crop).nitrogen.harvest_export,
        crop_prognostic(crop).nitrogen.storage,
        crop_prognostic(crop).nitrogen.leaf,
        crop_prognostic(crop).nitrogen.pool,
        crop_prognostic(crop).nitrogen.root,
        soil_carbon_fluxes(soil).input,
        soil_nitrogen_fluxes(soil).input,
        soil_water_prognostic(soil).storage,
        soil_properties(soil).layer_depth,
        residue_frac,
        output.crop.yield,
        output.calendar.harvest_date,
        output.calendar.harvesting_year,
        annual_row,
        day,
    )

    daily_sources = (
        crop = (
            growing_mask = crop_prognostic(crop).phenology.is_growing,
            storage_carbon = crop_prognostic(crop).carbon.storage,
        ),
        calendar = (
            harvesting_mask = crop_events(crop).harvest,
            sowing_event = crop_events(crop).sowing,
            harvest_event = crop_events(crop).harvest,
        ),
    )
    for (container_name, sources) in pairs(daily_sources)
        container = getproperty(output, container_name)
        for (field, source) in pairs(sources)
            if output_row === nothing
                setproperty!(
                    container,
                    field,
                    _append_output_row(getproperty(container, field), source),
                )
            else
                _write_output_row!(getproperty(container, field), output_row, source)
            end
        end
    end
    if day == 365 && annual_output_row === nothing
        annual_yield = max.(output.annual.yield, zero(T))
        harvesting_year = ifelse.(annual_yield .!= zero(T), Int32(1), Int32(0))
        output.calendar.harvest_date = _append_output_row(
            output.calendar.harvest_date, output.annual.harvest_date,
        )
        output.crop.yield = _append_output_row(output.crop.yield, annual_yield)
        output.calendar.harvesting_year = _append_output_row(
            output.calendar.harvesting_year, harvesting_year,
        )
        launch_1D!(
            reset_annual_harvest_kernel!,
            output.annual.yield,
            output.annual.harvest_date,
        )
    end
end

"""
    terminate_failed_crop!(crop, soil, output, residue_fraction, day; output_row, annual_output_row)

Harvest and remove crop stands flagged by the LPJmL negative-biomass test.
Normal calendar harvests have already set `is_growing` to zero and are ignored.
"""
function terminate_failed_crop!(crop,
                                soil,
                                output::Output,
                                residue_fraction::AbstractVector{T},
                                day::Integer;
                                output_row::Integer,
                                annual_output_row::Union{Nothing, Integer} = nothing,
) where {T <: AbstractFloat}
    day == 365 && annual_output_row === nothing &&
        throw(ArgumentError("annual_output_row is required on day 365"))
    annual_row = something(annual_output_row, 0)
    launch_1D!(
        terminate_failed_crop_kernel!,
        crop_events(crop).harvest,
        output.annual.harvest_date,
        crop_prognostic(crop).phenology.is_growing,
        crop_prognostic(crop).phenology.harvesting,
        crop_prognostic(crop).phenology.harvesting_previous,
        crop_fluxes(crop).carbon.yield,
        crop_fluxes(crop).carbon.harvest_export,
        output.annual.yield,
        crop_prognostic(crop).carbon.biomass,
        crop_prognostic(crop).carbon.storage,
        crop_prognostic(crop).carbon.leaf,
        crop_prognostic(crop).carbon.pool,
        crop_prognostic(crop).carbon.root,
        crop_fluxes(crop).nitrogen.harvest_export,
        crop_prognostic(crop).nitrogen.total,
        crop_prognostic(crop).nitrogen.storage,
        crop_prognostic(crop).nitrogen.leaf,
        crop_prognostic(crop).nitrogen.pool,
        crop_prognostic(crop).nitrogen.root,
        crop_prognostic(crop).nitrogen.pending_manure,
        crop_prognostic(crop).nitrogen.pending_fertilizer,
        crop_prognostic(crop).canopy.lai,
        crop_prognostic(crop).canopy.lai_npp_deficit,
        crop_canopy_auxiliary(crop).actual_lai,
        soil_carbon_fluxes(soil).input,
        soil_nitrogen_fluxes(soil).input,
        residue_fraction,
        output.crop.biomass,
        output.crop.lai,
        output.crop.storage_carbon,
        output.crop.growing_mask,
        output.calendar.harvesting_mask,
        output.calendar.harvest_event,
        output.crop.yield,
        output.calendar.harvest_date,
        output.calendar.harvesting_year,
        output_row,
        annual_row,
        day,
    )
    return nothing
end

@kernel inbounds = true function terminate_failed_crop_kernel!(
    harvest_event::AbstractVector{S},
    harvest_date::AbstractVector{S},
    is_growing::AbstractVector{S},
    harvesting::AbstractVector{B},
    harvesting_previous::AbstractVector{B},
    crop_yield::AbstractVector{T},
    carbon_harvest_export::AbstractVector{T},
    annual_yield::AbstractVector{T},
    biomass::AbstractVector{T},
    storage_carbon::AbstractVector{T},
    leaf_carbon::AbstractVector{T},
    pool_carbon::AbstractVector{T},
    root_carbon::AbstractVector{T},
    nitrogen_harvest_export::AbstractVector{T},
    total_nitrogen::AbstractVector{T},
    storage_nitrogen::AbstractVector{T},
    leaf_nitrogen::AbstractVector{T},
    pool_nitrogen::AbstractVector{T},
    root_nitrogen::AbstractVector{T},
    pending_manure::AbstractVector{T},
    pending_fertilizer::AbstractVector{T},
    lai::AbstractVector{T},
    lai_npp_deficit::AbstractVector{T},
    actual_lai::AbstractVector{T},
    carbon_input::AbstractMatrix{T},
    nitrogen_input::AbstractMatrix{T},
    residue_fraction::AbstractVector{T},
    output_biomass::AbstractMatrix{T},
    output_lai::AbstractMatrix{T},
    output_storage_carbon::AbstractMatrix{T},
    output_growing_mask::AbstractMatrix{S},
    output_harvesting_mask::AbstractMatrix{S},
    output_harvest_event::AbstractMatrix{S},
    output_yield::AbstractMatrix{T},
    output_harvest_date::AbstractMatrix{S},
    output_harvesting_year::AbstractMatrix{S},
    output_row::Integer,
    annual_output_row::Integer,
    day::Integer,
) where {T <: AbstractFloat, S <: Integer, B <: Bool}
    cell = @index(Global)
    failed = harvest_event[cell] != 0 && is_growing[cell] != 0
    if failed
        crop_yield[cell] = storage_carbon[cell]
        annual_yield[cell] += crop_yield[cell]
        harvest_date[cell] = S(day)
        carbon_harvest_export[cell] = crop_yield[cell] +
            (leaf_carbon[cell] + pool_carbon[cell]) *
            (one(T) - residue_fraction[cell])
        nitrogen_harvest_export[cell] = storage_nitrogen[cell] +
            (leaf_nitrogen[cell] + pool_nitrogen[cell]) *
            (one(T) - residue_fraction[cell])
        carbon_input[SURFACE_LITTER, cell] =
            (leaf_carbon[cell] + pool_carbon[cell]) * residue_fraction[cell]
        carbon_input[ROOT_LITTER, cell] = root_carbon[cell]
        nitrogen_input[SURFACE_LITTER, cell] =
            (leaf_nitrogen[cell] + pool_nitrogen[cell]) * residue_fraction[cell]
        nitrogen_input[ROOT_LITTER, cell] = root_nitrogen[cell]
        carbon_input[INCORPORATED_LITTER, cell] = zero(T)
        nitrogen_input[INCORPORATED_LITTER, cell] = zero(T)

        is_growing[cell] = zero(S)
        harvesting[cell] = false
        harvesting_previous[cell] = false
        biomass[cell] = zero(T)
        storage_carbon[cell] = zero(T)
        leaf_carbon[cell] = zero(T)
        pool_carbon[cell] = zero(T)
        root_carbon[cell] = zero(T)
        total_nitrogen[cell] = zero(T)
        storage_nitrogen[cell] = zero(T)
        leaf_nitrogen[cell] = zero(T)
        pool_nitrogen[cell] = zero(T)
        root_nitrogen[cell] = zero(T)
        pending_manure[cell] = zero(T)
        pending_fertilizer[cell] = zero(T)
        lai[cell] = zero(T)
        lai_npp_deficit[cell] = zero(T)
        actual_lai[cell] = zero(T)
    end
    size(output_biomass, 1) != 0 && (output_biomass[output_row, cell] = biomass[cell])
    size(output_lai, 1) != 0 && (output_lai[output_row, cell] = actual_lai[cell])
    size(output_storage_carbon, 1) != 0 &&
        (output_storage_carbon[output_row, cell] = storage_carbon[cell])
    size(output_growing_mask, 1) != 0 &&
        (output_growing_mask[output_row, cell] = is_growing[cell])
    if size(output_harvesting_mask, 1) != 0
        output_harvesting_mask[output_row, cell] = max(
            output_harvesting_mask[output_row, cell], failed ? one(S) : zero(S),
        )
    end
    if size(output_harvest_event, 1) != 0
        output_harvest_event[output_row, cell] = max(
            output_harvest_event[output_row, cell], failed ? one(S) : zero(S),
        )
    end

    # The regular harvest kernel has already emitted and cleared day 365.
    if day == 365
        if failed
            size(output_yield, 1) != 0 &&
                (output_yield[annual_output_row, cell] += crop_yield[cell])
            size(output_harvest_date, 1) != 0 &&
                (output_harvest_date[annual_output_row, cell] = S(day))
            size(output_harvesting_year, 1) != 0 &&
                (output_harvesting_year[annual_output_row, cell] = one(S))
        end
        annual_yield[cell] = zero(T)
        harvest_date[cell] = zero(S)
    end
    # Keep only failure events in the live state so the second residue-routing
    # pass cannot route a normal calendar harvest twice.
    harvest_event[cell] = failed ? one(S) : zero(S)
end

@kernel inbounds = true function harvest_state_kernel!(
    harvest_event::AbstractVector{S},
    harvest_date::AbstractVector{S},
    harvesting_previous::AbstractVector{B},
    harvesting::AbstractVector{B},
    is_growing::AbstractVector{S},
    crop_yield::AbstractVector{T},
    carbon_harvest_export::AbstractVector{T},
    annual_yield::AbstractVector{T},
    storage_carbon::AbstractVector{T},
    leaf_carbon::AbstractVector{T},
    pool_carbon::AbstractVector{T},
    root_carbon::AbstractVector{T},
    harvest_nitrogen::AbstractVector{T},
    storage_nitrogen::AbstractVector{T},
    leaf_nitrogen::AbstractVector{T},
    pool_nitrogen::AbstractVector{T},
    root_nitrogen::AbstractVector{T},
    carbon_input::AbstractMatrix{T},
    nitrogen_input::AbstractMatrix{T},
    soil_water_storage::AbstractMatrix{T},
    soil_layer_depth::AbstractVector{T},
    residue_fraction::AbstractVector{T},
    output_yield::AbstractMatrix{T},
    output_harvest_date::AbstractMatrix{S},
    output_harvesting_year::AbstractMatrix{S},
    annual_output_row::Integer,
    day::Integer,
) where {T <: AbstractFloat, S <: Integer, B <: Bool}
    cell = @index(Global)
    harvested = !harvesting_previous[cell] && harvesting[cell]
    event = harvested ? one(S) : zero(S)
    harvest_event[cell] = event
    if harvested
        harvest_date[cell] = S(day)
        is_growing[cell] = zero(S)
        crop_yield[cell] = storage_carbon[cell]
        annual_yield[cell] += crop_yield[cell]
        carbon_harvest_export[cell] = crop_yield[cell] +
            (leaf_carbon[cell] + pool_carbon[cell]) *
            (one(T) - residue_fraction[cell])
        harvest_nitrogen[cell] = storage_nitrogen[cell] +
            (leaf_nitrogen[cell] + pool_nitrogen[cell]) *
            (one(T) - residue_fraction[cell])
        carbon_input[SURFACE_LITTER, cell] =
            (leaf_carbon[cell] + pool_carbon[cell]) * residue_fraction[cell]
        carbon_input[ROOT_LITTER, cell] = root_carbon[cell]
        nitrogen_input[SURFACE_LITTER, cell] =
            (leaf_nitrogen[cell] + pool_nitrogen[cell]) * residue_fraction[cell]
        nitrogen_input[ROOT_LITTER, cell] = root_nitrogen[cell]
    else
        crop_yield[cell] = zero(T)
        carbon_harvest_export[cell] = zero(T)
        harvest_nitrogen[cell] = zero(T)
        carbon_input[SURFACE_LITTER, cell] = zero(T)
        carbon_input[ROOT_LITTER, cell] = zero(T)
        nitrogen_input[SURFACE_LITTER, cell] = zero(T)
        nitrogen_input[ROOT_LITTER, cell] = zero(T)
    end
    carbon_input[INCORPORATED_LITTER, cell] = zero(T)
    nitrogen_input[INCORPORATED_LITTER, cell] = zero(T)
    if day == 365 && annual_output_row != 0
        emitted_yield = max(annual_yield[cell], zero(T))
        size(output_yield, 1) != 0 &&
            (output_yield[annual_output_row, cell] = emitted_yield)
        size(output_harvest_date, 1) != 0 &&
            (output_harvest_date[annual_output_row, cell] = harvest_date[cell])
        size(output_harvesting_year, 1) != 0 &&
            (output_harvesting_year[annual_output_row, cell] =
                emitted_yield != zero(T) ? one(S) : zero(S))
        annual_yield[cell] = zero(T)
        harvest_date[cell] = zero(S)
    end
end

@kernel inbounds = true function reset_annual_harvest_kernel!(
    annual_yield::AbstractVector{T},
    harvest_date::AbstractVector{S},
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    annual_yield[cell] = zero(T)
    harvest_date[cell] = zero(S)
end
