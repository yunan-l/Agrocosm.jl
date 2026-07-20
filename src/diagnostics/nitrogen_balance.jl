"""
    NitrogenBalance

Optional daily nitrogen-budget diagnostics. All fields have shape
`(number_of_days, number_of_cells)` and use the model nitrogen unit (g N m⁻²).

The tracked system stock is plant total N plus soil mineral and organic N.
Organ N pools are derived partitions of plant total N and are not counted
again. Positive `residual` denotes nitrogen entering the tracked system
without appearing in final storage or a recorded boundary loss.
"""
mutable struct NitrogenBalance{M <: AbstractArray{<:AbstractFloat}}
    plant_before::M
    plant_after::M
    mineral_before::M
    mineral_after::M
    organic_before::M
    organic_after::M
    total_before::M
    total_after::M
    root_uptake::M
    seed_input::M
    prescribed_fertilizer_input::M
    prescribed_manure_input::M
    automatic_fertilizer_input::M
    harvest_export::M
    mineralization::M
    immobilization::M
    nitrification::M
    n2o_nitrification::M
    denitrification::M
    n2o_denitrification::M
    n2_denitrification::M
    volatilization::M
    gaseous_loss::M
    leaching_loss::M
    residual::M
    relative_residual::M
end

function init_nitrogen_balance(number_of_days::Integer,
                               number_of_cells::Integer,
                               device = identity;
                               T::Type{<:AbstractFloat} = Float32)
    allocate() = device(zeros(T, number_of_days, number_of_cells))
    return NitrogenBalance(
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(), allocate(), allocate(), allocate(), allocate(),
        allocate(),
    )
end

function record_nitrogen_balance_start!(balance::NitrogenBalance,
                                        day_index::Integer,
                                        crop::Crop,
                                        soil::Soil)
    @views begin
        balance.plant_before[day_index, :] .= crop.state.nitrogen.total
        balance.mineral_before[day_index, :] .= vec(sum(
            soil.nitrogen.nitrate .+ soil.nitrogen.ammonium; dims = 1,
        ))
        balance.organic_before[day_index, :] .= vec(sum(
            soil.nitrogen.litter; dims = 1,
        )) .+ vec(sum(soil.nitrogen.fast .+ soil.nitrogen.slow; dims = 1))
        balance.total_before[day_index, :] .=
            balance.plant_before[day_index, :] .+
            balance.mineral_before[day_index, :] .+
            balance.organic_before[day_index, :]
    end
    return nothing
end

function record_nitrogen_balance_end!(balance::NitrogenBalance,
                                      day_index::Integer,
                                      crop::Crop,
                                      soil::Soil)
    @views begin
        balance.plant_after[day_index, :] .= crop.state.nitrogen.total
        balance.mineral_after[day_index, :] .= vec(sum(
            soil.nitrogen.nitrate .+ soil.nitrogen.ammonium; dims = 1,
        ))
        balance.organic_after[day_index, :] .= vec(sum(
            soil.nitrogen.litter; dims = 1,
        )) .+ vec(sum(soil.nitrogen.fast .+ soil.nitrogen.slow; dims = 1))
        balance.total_after[day_index, :] .=
            balance.plant_after[day_index, :] .+
            balance.mineral_after[day_index, :] .+
            balance.organic_after[day_index, :]

        balance.root_uptake[day_index, :] .=
            crop.fluxes.nitrogen.uptake .- crop.fluxes.nitrogen.auto_fertilizer
        balance.seed_input[day_index, :] .= crop.fluxes.nitrogen.seed_input
        balance.prescribed_fertilizer_input[day_index, :] .=
            crop.fluxes.nitrogen.prescribed_fertilizer_input
        balance.prescribed_manure_input[day_index, :] .=
            crop.fluxes.nitrogen.prescribed_manure_input
        balance.automatic_fertilizer_input[day_index, :] .=
            crop.fluxes.nitrogen.auto_fertilizer
        balance.harvest_export[day_index, :] .= crop.fluxes.nitrogen.harvest_export
        balance.mineralization[day_index, :] .=
            vec(sum(soil.nitrogen.mineralization; dims = 1))
        balance.immobilization[day_index, :] .=
            vec(sum(soil.nitrogen.immobilization; dims = 1))
        balance.nitrification[day_index, :] .=
            vec(sum(soil.nitrogen.nitrification; dims = 1))
        balance.n2o_nitrification[day_index, :] .=
            vec(sum(soil.nitrogen.n2o_nitrification; dims = 1))
        balance.denitrification[day_index, :] .=
            vec(sum(soil.nitrogen.denitrification; dims = 1))
        balance.n2o_denitrification[day_index, :] .=
            vec(sum(soil.nitrogen.n2o_denitrification; dims = 1))
        balance.n2_denitrification[day_index, :] .=
            vec(sum(soil.nitrogen.n2_denitrification; dims = 1))
        balance.volatilization[day_index, :] .= soil.nitrogen.volatilization
        balance.gaseous_loss[day_index, :] .=
            balance.n2o_nitrification[day_index, :] .+
            balance.n2o_denitrification[day_index, :] .+
            balance.n2_denitrification[day_index, :] .+
            balance.volatilization[day_index, :]
        balance.leaching_loss[day_index, :] .= soil.nitrogen.leaching

        balance.residual[day_index, :] .=
            balance.total_before[day_index, :] .+
            balance.seed_input[day_index, :] .+
            balance.prescribed_fertilizer_input[day_index, :] .+
            balance.prescribed_manure_input[day_index, :] .+
            balance.automatic_fertilizer_input[day_index, :] .-
            balance.harvest_export[day_index, :] .-
            balance.gaseous_loss[day_index, :] .-
            balance.leaching_loss[day_index, :] .-
            balance.total_after[day_index, :]
        balance.relative_residual[day_index, :] .=
            balance.residual[day_index, :] ./ max.(
                abs.(balance.total_before[day_index, :]) .+
                balance.seed_input[day_index, :] .+
                balance.prescribed_fertilizer_input[day_index, :] .+
                balance.prescribed_manure_input[day_index, :] .+
                balance.automatic_fertilizer_input[day_index, :],
                eps(eltype(balance.residual)),
            )
    end
    return nothing
end
