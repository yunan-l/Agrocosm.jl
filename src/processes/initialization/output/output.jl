"""Time-series outputs for crop processes."""
mutable struct CropOutput{A, I, M}
    gpp::A                 # Daily gross primary production (gC m⁻² day⁻¹).
    npp::A                 # Daily net primary production after all plant respiration (gC m⁻² day⁻¹).
    lambda::A              # Optimal intercellular-to-ambient CO₂ ratio used by photosynthesis (0–1).
    potential_vcmax::A     # Water- and N-unlimited maximum carboxylation capacity (gC m⁻² day⁻¹).
    vcmax::A               # Realized maximum carboxylation capacity after limitations (gC m⁻² day⁻¹).
    nitrogen_limitation::A # Realized-to-potential `vcmax` ratio (0–1).
    respiration::A         # Total plant maintenance plus growth respiration (gC m⁻² day⁻¹).
    biomass::A             # Total live crop carbon stock (gC m⁻²).
    lai::A                 # Actual nonnegative leaf-area index (m² leaf m⁻² ground).
    storage_carbon::A      # Carbon in the harvestable storage organ (gC m⁻²).
    yield::A               # Annual harvested storage-organ carbon (gC m⁻² yr⁻¹).
    vegetation_carbon::M   # Daily leaf/root/pool/storage carbon stocks (gC m⁻²).
    vegetation_nitrogen::M # Daily leaf/root/pool/storage nitrogen contents (gN m⁻²).
    fphu::A                # Fraction of potential heat units accumulated (0–1+).
    water_deficit::A       # Daily crop water-deficit factor (0–100%).
    growing_mask::I        # Integer mask: one while a crop stand is active, otherwise zero.
end

"""Time-series outputs for soil stocks and fluxes."""
mutable struct SoilOutput{A, M}
    ecosystem_respiration::A   # Plant plus heterotrophic respiration (gC m⁻² day⁻¹).
    litter_carbon::M           # Carbon stocks in surface/incorporated/root litter (gC m⁻²).
    fast_carbon::M             # Fast soil-organic-carbon stock by layer (gC m⁻²).
    slow_carbon::M             # Slow soil-organic-carbon stock by layer (gC m⁻²).
    water_storage::M           # Liquid soil-water storage by layer (mm).
    litter_nitrogen::M         # Nitrogen stocks in the three litter classes (gN m⁻²).
    fast_nitrogen::M           # Fast soil-organic-nitrogen stock by layer (gN m⁻²).
    slow_nitrogen::M           # Slow soil-organic-nitrogen stock by layer (gN m⁻²).
    heterotrophic_respiration::A # Litter and soil respiration (gC m⁻² day⁻¹).
    evapotranspiration::A      # Soil evaporation plus crop transpiration (mm day⁻¹).
end

"""Time-series outputs for climate forcing and potential evaporation."""
mutable struct ClimateOutput{A}
    equilibrium_evapotranspiration::A # Priestley–Taylor equilibrium evaporation demand (mm day⁻¹).
    precipitation::A                 # Daily precipitation forcing (mm day⁻¹).
    temperature::A                   # Daily near-surface air temperature forcing (°C).
end

"""Time-series crop calendar diagnostics."""
mutable struct CalendarOutput{I}
    harvesting_mask::I # Daily mask indicating the harvest window/condition.
    harvesting_year::I # Simulation year associated with each annual harvest record.
    harvest_date::I    # Day of year of the recorded annual harvest.
    sowing_event::I    # Daily one-day sowing event indicator (0/1).
    harvest_event::I   # Daily one-day harvest event indicator (0/1).
end

"""In-progress annual crop outputs retained until the calendar-year boundary."""
mutable struct AnnualOutputAccumulator{A, I}
    yield::A        # Harvested storage carbon accumulated in the current output year (gC m⁻²).
    harvest_date::I # Latest harvest day in the current output year (1–365; 0 if absent).
end

"""Process-grouped model output container."""
mutable struct Output{C, S, F, K, A}
    crop::C     # Crop daily and annual output time series.
    soil::S     # Soil stock and flux output time series.
    climate::F  # Selected climate-forcing output time series.
    calendar::K # Sowing and harvest calendar output time series.
    annual::A   # In-progress annual output records required before year-end emission.
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
    annual = AnnualOutputAccumulator(
        device(zeros(T, cell_size)), device(zeros(Int32, cell_size)),
    )
    return Output(crop, soil, climate, calendar, annual)
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
