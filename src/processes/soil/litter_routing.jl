"""
    litter_tillage!(soil, crop_cal)

On sowing days, transfer the configured fraction of existing surface litter
to the incorporated litter pool. This follows LPJmL's `cultivate.c` ->
`tillage.c` order.
"""
function litter_tillage_reference!(soil::Soil, crop_cal::CropCalendar)
    callback = reshape(crop_cal.sowing_callback, (1, :))

    surface_carbon_before = copy(@view soil.carbon.litter[SURFACE_LITTER, :])
    surface_nitrogen_before = copy(@view soil.nitrogen.litter[SURFACE_LITTER, :])

    tilled_carbon = soil.management.tillage_fraction * soil.carbon.litter
    tilled_nitrogen = soil.management.tillage_fraction * soil.nitrogen.litter
    soil.carbon.litter .= soil.carbon.litter .* (1 .- callback) .+ tilled_carbon .* callback
    soil.nitrogen.litter .= soil.nitrogen.litter .* (1 .- callback) .+ tilled_nitrogen .* callback

    soil.management.tillage_carbon .=
        max.(surface_carbon_before .- @view(soil.carbon.litter[SURFACE_LITTER, :]), zero(eltype(surface_carbon_before)))
    soil.management.tillage_nitrogen .=
        max.(surface_nitrogen_before .- @view(soil.nitrogen.litter[SURFACE_LITTER, :]), zero(eltype(surface_nitrogen_before)))
    return nothing
end

function litter_tillage!(soil::Soil, crop_cal::CropCalendar)
    launch_custom!(
        litter_tillage_kernel!,
        soil.carbon.litter,
        size(soil.carbon.litter, 2),
        soil.nitrogen.litter,
        soil.management.tillage_fraction,
        crop_cal.sowing_callback,
        soil.management.tillage_carbon,
        soil.management.tillage_nitrogen,
    )
    return nothing
end

@kernel inbounds = true function litter_tillage_kernel!(
    carbon_litter::AbstractMatrix{T},
    nitrogen_litter::AbstractMatrix{T},
    tillage_fraction::AbstractMatrix{T},
    sowing_callback::AbstractVector{S},
    tillage_carbon::AbstractVector{T},
    tillage_nitrogen::AbstractVector{T},
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    callback = sowing_callback[cell]
    if callback != 0
        carbon_1 = carbon_litter[1, cell]
        carbon_2 = carbon_litter[2, cell]
        carbon_3 = carbon_litter[3, cell]
        nitrogen_1 = nitrogen_litter[1, cell]
        nitrogen_2 = nitrogen_litter[2, cell]
        nitrogen_3 = nitrogen_litter[3, cell]
        for destination in 1:3
            carbon_routed = tillage_fraction[destination, 1] * carbon_1 +
                tillage_fraction[destination, 2] * carbon_2 +
                tillage_fraction[destination, 3] * carbon_3
            nitrogen_routed = tillage_fraction[destination, 1] * nitrogen_1 +
                tillage_fraction[destination, 2] * nitrogen_2 +
                tillage_fraction[destination, 3] * nitrogen_3
            carbon_litter[destination, cell] = carbon_routed
            nitrogen_litter[destination, cell] = nitrogen_routed
        end
        tillage_carbon[cell] = max(carbon_1 - carbon_litter[1, cell], zero(T))
        tillage_nitrogen[cell] = max(nitrogen_1 - nitrogen_litter[1, cell], zero(T))
    else
        tillage_carbon[cell] = zero(T)
        tillage_nitrogen[cell] = zero(T)
    end
end

"""
    litter_bioturbation!(soil; lpjmlparams=lpjmlparams)

Apply LPJmL's daily bioturbation transfer from surface (`agtop`) to
incorporated (`agsub`) litter. Carbon and nitrogen are moved together and the
operation is conservative for each cell.
"""
function litter_bioturbation_reference!(soil::Soil;
                                        lpjmlparams::LPJmLParams = lpjmlparams)
    fraction = lpjmlparams.bioturbate

    soil.management.bioturbation_carbon .= @view(soil.carbon.litter[SURFACE_LITTER, :]) .* fraction
    soil.management.bioturbation_nitrogen .= @view(soil.nitrogen.litter[SURFACE_LITTER, :]) .* fraction

    @views soil.carbon.litter[INCORPORATED_LITTER, :] .+= soil.management.bioturbation_carbon
    @views soil.carbon.litter[SURFACE_LITTER, :] .-= soil.management.bioturbation_carbon
    @views soil.nitrogen.litter[INCORPORATED_LITTER, :] .+= soil.management.bioturbation_nitrogen
    @views soil.nitrogen.litter[SURFACE_LITTER, :] .-= soil.management.bioturbation_nitrogen
    return nothing
end

function litter_bioturbation!(soil::Soil;
                              lpjmlparams::LPJmLParams = lpjmlparams)
    launch_custom!(
        litter_bioturbation_kernel!,
        soil.carbon.litter,
        size(soil.carbon.litter, 2),
        soil.nitrogen.litter,
        soil.management.bioturbation_carbon,
        soil.management.bioturbation_nitrogen,
        eltype(soil.carbon.litter)(lpjmlparams.bioturbate),
    )
    return nothing
end

@kernel inbounds = true function litter_bioturbation_kernel!(
    carbon_litter::AbstractMatrix{T},
    nitrogen_litter::AbstractMatrix{T},
    bioturbation_carbon::AbstractVector{T},
    bioturbation_nitrogen::AbstractVector{T},
    fraction::T,
) where {T <: AbstractFloat}
    cell = @index(Global)
    carbon_flux = carbon_litter[SURFACE_LITTER, cell] * fraction
    nitrogen_flux = nitrogen_litter[SURFACE_LITTER, cell] * fraction
    bioturbation_carbon[cell] = carbon_flux
    bioturbation_nitrogen[cell] = nitrogen_flux
    carbon_litter[INCORPORATED_LITTER, cell] += carbon_flux
    carbon_litter[SURFACE_LITTER, cell] -= carbon_flux
    nitrogen_litter[INCORPORATED_LITTER, cell] += nitrogen_flux
    nitrogen_litter[SURFACE_LITTER, cell] -= nitrogen_flux
end

"""
Route today's harvested carbon residues through LPJmL's post-harvest tillage.

`harvest_crop.c` first adds shoot residues to `agtop` and roots to `bg`.
The harvested stand is then marked `KILL`; the same day's `killstand()` calls
`setaside()`, which calls `tillage()` when tillage is enabled. The root pool is
unchanged by the tillage matrix.
"""
function route_harvest_carbon_input_reference!(soil::Soil, crop_cal::CropCalendar)
    callback = reshape(crop_cal.harvest_callback, (1, :))
    litter_with_input = soil.carbon.litter .+
                        max.(soil.carbon.input, zero(eltype(soil.carbon.input)))
    routed_litter = soil.management.tillage_fraction * litter_with_input
    soil.management.tillage_carbon .+=
        max.((@view litter_with_input[SURFACE_LITTER, :]) .-
             (@view routed_litter[SURFACE_LITTER, :]),
             zero(eltype(litter_with_input))) .* vec(crop_cal.harvest_callback)
    soil.carbon.litter .= soil.carbon.litter .* (1 .- callback) .+
                          routed_litter .* callback
    return nothing
end

function route_harvest_carbon_input!(soil::Soil, crop_cal::CropCalendar)
    launch_custom!(
        route_harvest_litter_kernel!,
        soil.carbon.litter,
        size(soil.carbon.litter, 2),
        soil.carbon.input,
        soil.management.tillage_fraction,
        crop_cal.harvest_callback,
        soil.management.tillage_carbon,
    )
    return nothing
end

"""Route today's harvested nitrogen residues through post-harvest tillage."""
function route_harvest_nitrogen_input_reference!(soil::Soil, crop_cal::CropCalendar)
    callback = reshape(crop_cal.harvest_callback, (1, :))
    litter_with_input = soil.nitrogen.litter .+
                        max.(soil.nitrogen.input, zero(eltype(soil.nitrogen.input)))
    routed_litter = soil.management.tillage_fraction * litter_with_input
    soil.management.tillage_nitrogen .+=
        max.((@view litter_with_input[SURFACE_LITTER, :]) .-
             (@view routed_litter[SURFACE_LITTER, :]),
             zero(eltype(litter_with_input))) .* vec(crop_cal.harvest_callback)
    soil.nitrogen.litter .= soil.nitrogen.litter .* (1 .- callback) .+
                            routed_litter .* callback
    return nothing
end


function route_harvest_nitrogen_input!(soil::Soil, crop_cal::CropCalendar)
    launch_custom!(
        route_harvest_litter_kernel!,
        soil.nitrogen.litter,
        size(soil.nitrogen.litter, 2),
        soil.nitrogen.input,
        soil.management.tillage_fraction,
        crop_cal.harvest_callback,
        soil.management.tillage_nitrogen,
    )
    return nothing
end

@kernel inbounds = true function route_harvest_litter_kernel!(
    litter::AbstractMatrix{T},
    litter_input::AbstractMatrix{T},
    tillage_fraction::AbstractMatrix{T},
    harvest_callback::AbstractVector{S},
    tillage_flux::AbstractVector{T},
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    if harvest_callback[cell] != 0
        litter_1 = litter[1, cell] + max(litter_input[1, cell], zero(T))
        litter_2 = litter[2, cell] + max(litter_input[2, cell], zero(T))
        litter_3 = litter[3, cell] + max(litter_input[3, cell], zero(T))
        routed_1 = tillage_fraction[1, 1] * litter_1 +
            tillage_fraction[1, 2] * litter_2 + tillage_fraction[1, 3] * litter_3
        tillage_flux[cell] += max(litter_1 - routed_1, zero(T))
        litter[1, cell] = routed_1
        litter[2, cell] = tillage_fraction[2, 1] * litter_1 +
            tillage_fraction[2, 2] * litter_2 + tillage_fraction[2, 3] * litter_3
        litter[3, cell] = tillage_fraction[3, 1] * litter_1 +
            tillage_fraction[3, 2] * litter_2 + tillage_fraction[3, 3] * litter_3
    end
end
