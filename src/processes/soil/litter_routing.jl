"""
    litter_tillage!(soil, crop_cal)

On sowing days, transfer the configured fraction of existing surface litter
to the incorporated litter pool. This follows LPJmL's `cultivate.c` ->
`tillage.c` order.
"""
function litter_tillage!(soil::Soil, crop_cal::CropCalendar)
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

"""
    litter_bioturbation!(soil; lpjmlparams=lpjmlparams)

Apply LPJmL's daily bioturbation transfer from surface (`agtop`) to
incorporated (`agsub`) litter. Carbon and nitrogen are moved together and the
operation is conservative for each cell.
"""
function litter_bioturbation!(soil::Soil;
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

"""
Route today's harvested carbon residues through LPJmL's post-harvest tillage.

`harvest_crop.c` first adds shoot residues to `agtop` and roots to `bg`.
The harvested stand is then marked `KILL`; the same day's `killstand()` calls
`setaside()`, which calls `tillage()` when tillage is enabled. The root pool is
unchanged by the tillage matrix.
"""
function route_harvest_carbon_input!(soil::Soil, crop_cal::CropCalendar)
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

"""Route today's harvested nitrogen residues through post-harvest tillage."""
function route_harvest_nitrogen_input!(soil::Soil, crop_cal::CropCalendar)
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
