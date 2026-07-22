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
        day,
    )

    # soil_carbon_fluxes(soil).input .= vcat(reshape((crop_prognostic(crop).carbon.leaf .+ crop_prognostic(crop).carbon.pool) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop_prognostic(crop).carbon.root, (1, :))) .* reshape(crop_events(crop).harvest, (1, :))
    # soil_nitrogen_fluxes(soil).input .= vcat(reshape((crop_prognostic(crop).nitrogen.leaf .+ crop_prognostic(crop).nitrogen.pool) .* residue_frac, (1, :)), device(zeros(Float32, (1, cell_size))), reshape(crop_prognostic(crop).nitrogen.root, (1, :))) .* reshape(crop_events(crop).harvest, (1, :))
    # idx = ((crop_prognostic(crop).phenology.harvesting_previous .== true) .& (crop_prognostic(crop).phenology.harvesting .== true)) .| ((crop_prognostic(crop).phenology.harvesting_previous .== true) .& (crop_prognostic(crop).phenology.harvesting .== false)) .| ((crop_prognostic(crop).phenology.harvesting_previous .== false) .& (crop_prognostic(crop).phenology.harvesting .== false))
    # crop_events(crop).harvest[idx] .= 0
    # crop_events(crop).harvest .= ifelse.(((crop_prognostic(crop).phenology.harvesting_previous .== true) .& (crop_prognostic(crop).phenology.harvesting .== true)) .| ((crop_prognostic(crop).phenology.harvesting_previous .== true) .& (crop_prognostic(crop).phenology.harvesting .== false)) .| ((crop_prognostic(crop).phenology.harvesting_previous .== false) .& (crop_prognostic(crop).phenology.harvesting .== false)), 0, crop_events(crop).harvest)

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
    if day == 365
        annual_yield = max.(output.annual.yield, zero(T))
        harvesting_year = ifelse.(annual_yield .!= zero(T), Int32(1), Int32(0))
        if annual_output_row === nothing
            output.calendar.harvest_date = _append_output_row(
                output.calendar.harvest_date, output.annual.harvest_date,
            )
            output.crop.yield = _append_output_row(
                output.crop.yield, annual_yield,
            )
            output.calendar.harvesting_year = _append_output_row(
                output.calendar.harvesting_year, harvesting_year,
            )
        else
            _write_output_row!(
                output.calendar.harvest_date, annual_output_row, output.annual.harvest_date,
            )
            _write_output_row!(output.crop.yield, annual_output_row, annual_yield)
            _write_output_row!(
                output.calendar.harvesting_year,
                annual_output_row,
                harvesting_year,
            )
        end
        output.annual.yield .= zero(T)
        output.annual.harvest_date .= 0
    end
    # crop_prognostic(crop).carbon.root = crop_prognostic(crop).carbon.root .* (1 .- crop_events(crop).harvest)
    # crop_prognostic(crop).carbon.leaf = crop_prognostic(crop).carbon.leaf .* (1 .- crop_events(crop).harvest)
    # crop_prognostic(crop).carbon.storage = crop_prognostic(crop).carbon.storage .* (1 .- crop_events(crop).harvest)
    # crop_prognostic(crop).carbon.pool = crop_prognostic(crop).carbon.pool .* (1 .- crop_events(crop).harvest)

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
end
