"""
harvest_crop!(crop_cal, crop, soil, day)

Handle harvest-day biomass removal, residue transfer, and crop state reset.
"""
function harvest_crop!(crop_cal::CropCalendar,
                       crop::Crop,
                       soil::Soil,
                       output::Output,
                       residue_frac::AbstractArray{T},
                       day::Int;
                       output_row::Union{Nothing, Integer} = nothing,
                       annual_output_row::Union{Nothing, Integer} = nothing
) where {T <: AbstractFloat}

    launch_1D!(
        harvest_state_kernel!,
        crop_cal.harvest_callback,
        crop_cal.harvest_date,
        crop.phenology.harvesting_previous,
        crop.phenology.harvesting,
        crop.phenology.is_growing,
        crop.carbon.yield,
        crop.carbon.storage,
        crop.carbon.leaf,
        crop.carbon.pool,
        crop.carbon.root,
        crop.nitrogen.harvest_export,
        crop.nitrogen.storage,
        crop.nitrogen.leaf,
        crop.nitrogen.pool,
        crop.nitrogen.root,
        soil.carbon.input,
        soil.nitrogen.input,
        residue_frac,
        day,
    )

    # Do not clear plant N here. harvest_crop! has already set is_growing = 0,
    # so the later same-day crop_nitrogen! call clears total N in
    # nuptake_crop! and all organ N pools in allocate_crop_nitrogen!. Keeping
    # that operation in the kernels avoids five redundant GPU broadcasts.
    # crop.nitrogen.total .*= one(T) .- crop_cal.harvest_callback
    # crop.nitrogen.leaf .*= one(T) .- crop_cal.harvest_callback
    # crop.nitrogen.root .*= one(T) .- crop_cal.harvest_callback
    # crop.nitrogen.storage .*= one(T) .- crop_cal.harvest_callback
    # crop.nitrogen.pool .*= one(T) .- crop_cal.harvest_callback

    # soil.carbon.input .= vcat(reshape((crop.carbon.leaf .+ crop.carbon.pool) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop.carbon.root, (1, :))) .* reshape(crop_cal.harvest_callback, (1, :))
    # soil.nitrogen.input .= vcat(reshape((crop.nitrogen.leaf .+ crop.nitrogen.pool) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop.nitrogen.root, (1, :))) .* reshape(crop_cal.harvest_callback, (1, :))
    # idx = ((crop.phenology.harvesting_previous .== true) .& (crop.phenology.harvesting .== true)) .| ((crop.phenology.harvesting_previous .== true) .& (crop.phenology.harvesting .== false)) .| ((crop.phenology.harvesting_previous .== false) .& (crop.phenology.harvesting .== false))
    # crop_cal.harvest_callback[idx] .= 0
    # crop_cal.harvest_callback .= ifelse.(((crop.phenology.harvesting_previous .== true) .& (crop.phenology.harvesting .== true)) .| ((crop.phenology.harvesting_previous .== true) .& (crop.phenology.harvesting .== false)) .| ((crop.phenology.harvesting_previous .== false) .& (crop.phenology.harvesting .== false)), 0, crop_cal.harvest_callback)

    # update harvesting variables
    launch_1D!(
        root_zone_mean_water_kernel!,
        crop.water.root_zone_water,
        soil.water.storage,
        soil.properties.layer_depth,
        3,
    )
    daily_sources = (
        crop = (
            growing_mask = crop.phenology.is_growing,
            storage_carbon = crop.carbon.storage,
            water_deficit = crop.water.root_zone_water,
        ),
        calendar = (
            harvesting_mask = crop_cal.harvest_callback,
            sowing_callback = crop_cal.sowing_callback,
            harvest_callback = crop_cal.harvest_callback,
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
        crop_cal.harvesting_year .= ifelse.(crop.carbon.yield .!= 0.0f0, 1, 0)
        crop.carbon.yield .= max.(crop.carbon.yield, 0.0f0)
        if annual_output_row === nothing
            output.calendar.harvest_date = _append_output_row(
                output.calendar.harvest_date, crop_cal.harvest_date,
            )
            output.crop.yield = _append_output_row(
                output.crop.yield, crop.carbon.yield,
            )
            output.calendar.harvesting_year = _append_output_row(
                output.calendar.harvesting_year, crop_cal.harvesting_year,
            )
        else
            _write_output_row!(
                output.calendar.harvest_date, annual_output_row, crop_cal.harvest_date,
            )
            _write_output_row!(output.crop.yield, annual_output_row, crop.carbon.yield)
            _write_output_row!(
                output.calendar.harvesting_year,
                annual_output_row,
                crop_cal.harvesting_year,
            )
        end
        crop.carbon.yield .= 0.0f0
        crop_cal.harvest_date .= 0
    end
    # crop.carbon.organs = crop.carbon.organs .* (1 .- reshape(crop_cal.harvest_callback, (1, :)))
    # crop.carbon.root = crop.carbon.root .* (1 .- crop_cal.harvest_callback)
    # crop.carbon.leaf = crop.carbon.leaf .* (1 .- crop_cal.harvest_callback)
    # crop.carbon.storage = crop.carbon.storage .* (1 .- crop_cal.harvest_callback)
    # crop.carbon.pool = crop.carbon.pool .* (1 .- crop_cal.harvest_callback)

end

@kernel inbounds = true function harvest_state_kernel!(
    harvest_callback::AbstractVector{S},
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
    callback = harvested ? one(S) : zero(S)
    harvest_callback[cell] = callback
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
