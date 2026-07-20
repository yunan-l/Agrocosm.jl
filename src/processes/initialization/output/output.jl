"""Time-series outputs for crop processes."""
mutable struct CropOutput{A, I, M}
    gpp::A
    npp::A
    lambda::A
    potential_vcmax::A
    vcmax::A
    nitrogen_limitation::A
    respiration::A
    biomass::A
    lai::A
    storage_carbon::A
    yield::A
    vegetation_carbon::M
    vegetation_nitrogen::M
    fphu::A
    water_deficit::A
    growing_mask::I
end

"""Time-series outputs for soil stocks and fluxes."""
mutable struct SoilOutput{A, M}
    ecosystem_respiration::A
    litter_carbon::M
    fast_carbon::M
    slow_carbon::M
    water_storage::M
    litter_nitrogen::M
    fast_nitrogen::M
    slow_nitrogen::M
    heterotrophic_respiration::A
    evapotranspiration::A
end

"""Time-series outputs for climate forcing and potential evaporation."""
mutable struct ClimateOutput{A}
    equilibrium_evapotranspiration::A
    precipitation::A
    temperature::A
end

"""Time-series crop calendar diagnostics."""
mutable struct CalendarOutput{I}
    harvesting_mask::I
    harvesting_year::I
    harvest_date::I
    sowing_event::I
    harvest_event::I
end

"""Process-grouped model output container."""
mutable struct Output{C, S, F, K}
    crop::C
    soil::S
    climate::F
    calendar::K
end

init_output(cell_size::Int, device; kwargs...) =
    init_output(Float32, cell_size, device; kwargs...)
function init_output(::Type{T},
                     cell_size::Int,
                     device;
                     vegc_pools::Int = 4,
                     litc_layers::Int = 3,
                     soil_layers::Int = 5) where {T <: AbstractFloat}
    # Output rows represent completed simulation steps only. Initial model
    # state lives in `crop`/`soil`; it is not a synthetic day-zero output.
    scalar_output() = device(zeros(T, 0, cell_size))
    integer_output() = device(zeros(Int32, 0, cell_size))

    crop = CropOutput(
        scalar_output(), scalar_output(), scalar_output(), scalar_output(),
        scalar_output(), scalar_output(), scalar_output(), scalar_output(),
        scalar_output(), scalar_output(), scalar_output(),
        device(zeros(T, 0, vegc_pools * cell_size)),
        device(zeros(T, 0, vegc_pools * cell_size)),
        scalar_output(), scalar_output(), integer_output(),
    )

    soil = SoilOutput(
        scalar_output(),
        device(zeros(T, 0, litc_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        device(zeros(T, 0, litc_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        device(zeros(T, 0, soil_layers * cell_size)),
        scalar_output(), scalar_output(),
    )

    climate = ClimateOutput(scalar_output(), scalar_output(), scalar_output())
    calendar = CalendarOutput(
        integer_output(), integer_output(), integer_output(),
        integer_output(), integer_output(),
    )
    return Output(crop, soil, climate, calendar)
end

"""Grow a backend array once for a simulation block, preserving existing rows."""
function _extend_output_rows(array::AbstractMatrix, additional_rows::Integer)
    additional_rows <= 0 && return array
    old_rows, columns = size(array)
    extended = similar(array, old_rows + additional_rows, columns)
    fill!(extended, zero(eltype(extended)))
    @views extended[1:old_rows, :] .= array
    return extended
end

"""
    prepare_output_block!(output, daily_rows, annual_rows)

Reserve all crop/calendar output rows once before a daily simulation block.
The returned indices point to the first newly allocated daily and annual rows.
"""
function prepare_output_block!(output::Output,
                               daily_rows::Integer,
                               annual_rows::Integer)
    first_daily_row = size(output.crop.gpp, 1) + 1
    first_annual_row = size(output.crop.yield, 1) + 1

    for field in (
        :gpp, :npp, :lambda, :potential_vcmax, :vcmax,
        :nitrogen_limitation, :respiration, :biomass, :lai,
        :storage_carbon, :fphu, :water_deficit, :growing_mask,
    )
        setproperty!(
            output.crop,
            field,
            _extend_output_rows(getproperty(output.crop, field), daily_rows),
        )
    end
    for field in (:harvesting_mask, :sowing_event, :harvest_event)
        setproperty!(
            output.calendar,
            field,
            _extend_output_rows(getproperty(output.calendar, field), daily_rows),
        )
    end

    output.crop.yield = _extend_output_rows(output.crop.yield, annual_rows)
    output.calendar.harvest_date =
        _extend_output_rows(output.calendar.harvest_date, annual_rows)
    output.calendar.harvesting_year =
        _extend_output_rows(output.calendar.harvesting_year, annual_rows)

    return (; first_daily_row, first_annual_row)
end

@inline function _write_output_row!(destination::AbstractMatrix,
                                    row::Integer,
                                    source::AbstractVector)
    @views destination[row:row, :] .= reshape(source, 1, :)
    return nothing
end

function _append_output_row(array::AbstractMatrix, source::AbstractVector)
    row = size(array, 1) + 1
    extended = _extend_output_rows(array, 1)
    _write_output_row!(extended, row, source)
    return extended
end
