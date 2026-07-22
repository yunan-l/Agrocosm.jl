"""
ClimateDataLoader(climate, data_index, device)

Extract climate forcing slices for selected grid points and years.
"""
function ClimateDataLoader(climate::NamedTuple, 
                           data_index::Vector{Int},
                           device;
                           T::Type{<:AbstractFloat} = Float32,
)

    loaded_climate = (
        temp_spinup = T.(climate.temp_spinup[:, data_index]),
        temp = T.(climate.temp[:, data_index]),
        prec = T.(climate.prec[:, data_index]),
        sw = T.(climate.swdown[:, data_index]),
        lw = T.(climate.lwnet[:, data_index]),
        co2 = T.(climate.co2),
    )

    # Normalize optional wind forcing to one internal field name while keeping
    # existing climate archives (which do not contain wind) fully compatible.
    if hasproperty(climate, :windspeed)
        loaded_climate = merge(
            loaded_climate,
            (wind = T.(climate.windspeed[:, data_index]),),
        )
    end

    if hasproperty(climate, :co2_daily)
        loaded_climate = merge(loaded_climate, (co2_daily = climate.co2_daily,))
    end

    return device(loaded_climate)
end
