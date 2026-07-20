"""
harvest_crop!(crop, soil, output, residue_fraction, day)

Handle harvest-day biomass removal, residue transfer, and crop state reset.
"""
function harvest_crop!(crop::Crop,
                       soil::Soil,
                       output::Output,
                       residue_frac::AbstractArray{T},
                       day::Int;
                       output_row::Union{Nothing, Integer} = nothing,
                       annual_output_row::Union{Nothing, Integer} = nothing
) where {T <: AbstractFloat}

    launch_1D!(
        harvest_state_kernel!,
        crop.events.harvest,
        crop.state.calendar.harvest_date,
        crop.state.phenology.harvesting_previous,
        crop.state.phenology.harvesting,
        crop.state.phenology.is_growing,
        crop.fluxes.carbon.yield,
        crop.state.carbon.storage,
        crop.state.carbon.leaf,
        crop.state.carbon.pool,
        crop.state.carbon.root,
        crop.fluxes.nitrogen.harvest_export,
        crop.state.nitrogen.storage,
        crop.state.nitrogen.leaf,
        crop.state.nitrogen.pool,
        crop.state.nitrogen.root,
        soil.carbon.input,
        soil.nitrogen.input,
        residue_frac,
        day,
    )

    # Do not clear plant N here. harvest_crop! has already set is_growing = 0,
    # so the later same-day crop_nitrogen! call clears total N in
    # nuptake_crop! and all organ N pools in allocate_crop_nitrogen!. Keeping
    # that operation in the kernels avoids five redundant GPU broadcasts.
    # crop.state.nitrogen.total .*= one(T) .- crop.events.harvest
    # crop.state.nitrogen.leaf .*= one(T) .- crop.events.harvest
    # crop.state.nitrogen.root .*= one(T) .- crop.events.harvest
    # crop.state.nitrogen.storage .*= one(T) .- crop.events.harvest
    # crop.state.nitrogen.pool .*= one(T) .- crop.events.harvest

    # soil.carbon.input .= vcat(reshape((crop.state.carbon.leaf .+ crop.state.carbon.pool) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop.state.carbon.root, (1, :))) .* reshape(crop.events.harvest, (1, :))
    # soil.nitrogen.input .= vcat(reshape((crop.state.nitrogen.leaf .+ crop.state.nitrogen.pool) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop.state.nitrogen.root, (1, :))) .* reshape(crop.events.harvest, (1, :))
    # idx = ((crop.state.phenology.harvesting_previous .== true) .& (crop.state.phenology.harvesting .== true)) .| ((crop.state.phenology.harvesting_previous .== true) .& (crop.state.phenology.harvesting .== false)) .| ((crop.state.phenology.harvesting_previous .== false) .& (crop.state.phenology.harvesting .== false))
    # crop.events.harvest[idx] .= 0
    # crop.events.harvest .= ifelse.(((crop.state.phenology.harvesting_previous .== true) .& (crop.state.phenology.harvesting .== true)) .| ((crop.state.phenology.harvesting_previous .== true) .& (crop.state.phenology.harvesting .== false)) .| ((crop.state.phenology.harvesting_previous .== false) .& (crop.state.phenology.harvesting .== false)), 0, crop.events.harvest)

    # update harvesting variables
    launch_1D!(
        root_zone_mean_water_kernel!,
        crop.auxiliary.stress.root_zone_water,
        soil.water.storage,
        soil.properties.layer_depth,
        3,
    )
    daily_sources = (
        crop = (
            growing_mask = crop.state.phenology.is_growing,
            storage_carbon = crop.state.carbon.storage,
            water_deficit = crop.auxiliary.stress.root_zone_water,
        ),
        calendar = (
            harvesting_mask = crop.events.harvest,
            sowing_event = crop.events.sowing,
            harvest_event = crop.events.harvest,
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
    if day == 365
        crop.state.calendar.harvesting_year .= ifelse.(crop.fluxes.carbon.yield .!= 0.0f0, 1, 0)
        crop.fluxes.carbon.yield .= max.(crop.fluxes.carbon.yield, 0.0f0)
        if annual_output_row === nothing
            output.calendar.harvest_date = _append_output_row(
                output.calendar.harvest_date, crop.state.calendar.harvest_date,
            )
            output.crop.yield = _append_output_row(
                output.crop.yield, crop.fluxes.carbon.yield,
            )
            output.calendar.harvesting_year = _append_output_row(
                output.calendar.harvesting_year, crop.state.calendar.harvesting_year,
            )
        else
            _write_output_row!(
                output.calendar.harvest_date, annual_output_row, crop.state.calendar.harvest_date,
            )
            _write_output_row!(output.crop.yield, annual_output_row, crop.fluxes.carbon.yield)
            _write_output_row!(
                output.calendar.harvesting_year,
                annual_output_row,
                crop.state.calendar.harvesting_year,
            )
        end
        crop.fluxes.carbon.yield .= 0.0f0
        crop.state.calendar.harvest_date .= 0
    end
    # crop.state.carbon.root = crop.state.carbon.root .* (1 .- crop.events.harvest)
    # crop.state.carbon.leaf = crop.state.carbon.leaf .* (1 .- crop.events.harvest)
    # crop.state.carbon.storage = crop.state.carbon.storage .* (1 .- crop.events.harvest)
    # crop.state.carbon.pool = crop.state.carbon.pool .* (1 .- crop.events.harvest)

end

@kernel inbounds = true function harvest_state_kernel!(
    harvest_event::AbstractVector{S},
    harvest_date::AbstractVector{S},
    harvesting_previous::AbstractVector{B},
    harvesting::AbstractVector{B},
    is_growing::AbstractVector{S},
    crop_yield::AbstractVector{T},
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
    residue_fraction::AbstractVector{T},
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
        harvest_nitrogen[cell] = zero(T)
        carbon_input[SURFACE_LITTER, cell] = zero(T)
        carbon_input[ROOT_LITTER, cell] = zero(T)
        nitrogen_input[SURFACE_LITTER, cell] = zero(T)
        nitrogen_input[ROOT_LITTER, cell] = zero(T)
    end
    carbon_input[INCORPORATED_LITTER, cell] = zero(T)
    nitrogen_input[INCORPORATED_LITTER, cell] = zero(T)
end

@kernel inbounds = true function root_zone_mean_water_kernel!(
    destination::AbstractVector{T},
    storage::AbstractMatrix{T},
    layer_depth::AbstractVector{T},
    root_layers::Integer,
) where {T <: AbstractFloat}
    cell = @index(Global)
    total = zero(T)
    for layer in 1:root_layers
        total += storage[layer, cell] / layer_depth[layer]
    end
    destination[cell] = total / T(root_layers)
end
