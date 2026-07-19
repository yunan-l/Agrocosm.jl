"""
ClimateDataLoader(climate, data_index, device)

Extract climate forcing slices for selected grid points and years.
"""
function ClimateDataLoader(climate::NamedTuple, 
                           data_index::Vector{Int},
                           device
)

    loaded_climate = (
        temp_spinup = climate.temp_spinup[:, data_index],
        temp = climate.temp[:, data_index],
        prec = climate.prec[:, data_index],
        sw = climate.swdown[:, data_index],
        lw = climate.lwnet[:, data_index],
        co2 = climate.co2
    )

    # Normalize optional wind forcing to one internal field name while keeping
    # existing climate archives (which do not contain wind) fully compatible.
    if hasproperty(climate, :windspeed)
        loaded_climate = merge(loaded_climate, (wind = climate.windspeed[:, data_index],))
    end

    return device(loaded_climate)
end
