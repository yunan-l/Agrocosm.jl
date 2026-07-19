"""
    CarbonBalance

Optional daily carbon-budget diagnostics. All fields have shape
`(number_of_days, number_of_cells)` and use the model carbon unit (g C m⁻²).

The tracked stock is crop organ carbon plus litter, fast-soil, and slow-soil
carbon. NPP, seed carbon, and manure carbon are boundary inputs; harvest
export and heterotrophic respiration are boundary outputs. Residue transfer
is recorded for inspection but is internal to the tracked system.
"""
mutable struct CarbonBalance{M <: AbstractArray{<:AbstractFloat}}
    plant_before::M
    plant_after::M
    soil_before::M
    soil_after::M
    total_before::M
    total_after::M
    net_primary_production::M
    seed_input::M
    manure_input::M
    residue_transfer::M
    harvest_export::M
    heterotrophic_respiration::M
    residual::M
    relative_residual::M
end

function init_carbon_balance(number_of_days::Integer,
                             number_of_cells::Integer,
                             device = identity;
                             T::Type{<:AbstractFloat} = Float32)
    allocate() = device(zeros(T, number_of_days, number_of_cells))
    return CarbonBalance(
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(),
    )
end

function crop_carbon_stock(crop::Crop)
    return crop.carbon.leaf .+ crop.carbon.root .+
           crop.carbon.pool .+ crop.carbon.storage
end

function soil_carbon_stock(soil::Soil)
    return vec(sum(soil.carbon.litter; dims = 1)) .+
           vec(sum(soil.carbon.fast .+ soil.carbon.slow; dims = 1))
end

function record_carbon_balance_start!(balance::CarbonBalance,
                                      day_index::Integer,
                                      crop::Crop,
                                      soil::Soil)
    @views begin
        balance.plant_before[day_index, :] .= crop_carbon_stock(crop)
        balance.soil_before[day_index, :] .= soil_carbon_stock(soil)
        balance.total_before[day_index, :] .=
            balance.plant_before[day_index, :] .+
            balance.soil_before[day_index, :]
    end
    return nothing
end

function record_carbon_balance_after_cultivate!(balance::CarbonBalance,
                                                day_index::Integer,
                                                crop::Crop;
                                                lpjmlparams::LPJmLParams = lpjmlparams)
    @views begin
        balance.seed_input[day_index, :] .=
            max.(crop_carbon_stock(crop) .-
                 balance.plant_before[day_index, :], zero(eltype(balance.seed_input))) .*
            crop.calendar.sowing_callback
        balance.manure_input[day_index, :] .=
            crop.nitrogen.prescribed_manure_input .* lpjmlparams.manure_cn
    end
    return nothing
end

function record_carbon_balance_after_harvest!(balance::CarbonBalance,
                                              day_index::Integer,
                                              crop::Crop,
                                              soil::Soil,
                                              residue_fraction)
    callback = crop.calendar.harvest_callback
    @views begin
        balance.residue_transfer[day_index, :] .=
            vec(sum(soil.carbon.input; dims = 1))
        balance.harvest_export[day_index, :] .=
            (crop.carbon.storage .+
             (crop.carbon.leaf .+ crop.carbon.pool) .* (one(eltype(residue_fraction)) .- residue_fraction)) .*
            callback
    end
    return nothing
end

function record_carbon_balance_end!(balance::CarbonBalance,
                                    day_index::Integer,
                                    crop::Crop,
                                    soil::Soil)
    @views begin
        balance.plant_after[day_index, :] .= crop_carbon_stock(crop)
        balance.soil_after[day_index, :] .= soil_carbon_stock(soil)
        balance.total_after[day_index, :] .=
            balance.plant_after[day_index, :] .+
            balance.soil_after[day_index, :]
        balance.net_primary_production[day_index, :] .= crop.carbon.npp
        balance.heterotrophic_respiration[day_index, :] .=
            soil.carbon.heterotrophic_respiration

        balance.residual[day_index, :] .=
            balance.total_before[day_index, :] .+
            balance.seed_input[day_index, :] .+
            balance.manure_input[day_index, :] .+
            balance.net_primary_production[day_index, :] .-
            balance.harvest_export[day_index, :] .-
            balance.heterotrophic_respiration[day_index, :] .-
            balance.total_after[day_index, :]
        balance.relative_residual[day_index, :] .=
            balance.residual[day_index, :] ./ max.(
                abs.(balance.total_before[day_index, :]) .+
                balance.seed_input[day_index, :] .+
                balance.manure_input[day_index, :] .+
                abs.(balance.net_primary_production[day_index, :]),
                eps(eltype(balance.residual)),
            )
    end
    return nothing
end
