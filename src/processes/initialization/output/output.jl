"""Time-series outputs for crop processes."""
mutable struct CropOutput{A, I, M}
    gpp::A
    npp::A
    lambda::A
    potential_vmax::A
    vmax::A
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
    sowing_callback::I
    harvest_callback::I
end

"""Process-grouped model output container."""
mutable struct Output{C, S, F, K}
    crop::C
    soil::S
    climate::F
    calendar::K
end

function init_output(cell_size::Int,
                     device;
                     vegc_pools::Int = 4,
                     litc_layers::Int = 3,
                     soil_layers::Int = 5)
    scalar_output() = device(zeros(Float32, 1, cell_size))
    integer_output() = device(zeros(Int32, 1, cell_size))

    crop = CropOutput(
        scalar_output(), scalar_output(), scalar_output(), scalar_output(),
        scalar_output(), scalar_output(), scalar_output(), scalar_output(),
        scalar_output(), scalar_output(), scalar_output(),
        device(zeros(Float32, 1, vegc_pools * cell_size)),
        device(zeros(Float32, 1, vegc_pools * cell_size)),
        scalar_output(), scalar_output(), integer_output(),
    )

    soil = SoilOutput(
        scalar_output(),
        device(zeros(Float32, 1, litc_layers * cell_size)),
        device(zeros(Float32, 1, soil_layers * cell_size)),
        device(zeros(Float32, 1, soil_layers * cell_size)),
        device(zeros(Float32, 1, soil_layers * cell_size)),
        device(zeros(Float32, 1, litc_layers * cell_size)),
        device(zeros(Float32, 1, soil_layers * cell_size)),
        device(zeros(Float32, 1, soil_layers * cell_size)),
        scalar_output(), scalar_output(),
    )

    climate = ClimateOutput(scalar_output(), scalar_output(), scalar_output())
    calendar = CalendarOutput(
        integer_output(), integer_output(), integer_output(),
        integer_output(), integer_output(),
    )
    return Output(crop, soil, climate, calendar)
end
