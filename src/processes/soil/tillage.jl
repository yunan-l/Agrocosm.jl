"""Apply the sowing-day reduction in topsoil bulk density caused by tillage."""
function tillage_hydraulics_reference!(soil, crop;
                                       lpjmlparams::LPJmLParams = lpjmlparams)
    event = reshape(crop_events(crop).sowing .!= 0, (1, :))
    density_factor = soil_management_prognostic(soil).tillage_density_factor
    tilled_factor = density_factor .-
        (density_factor .- eltype(density_factor)(0.667)) .* lpjmlparams.mixing_efficiency
    density_factor .= ifelse.(event, tilled_factor, density_factor)
    return nothing
end

function tillage_hydraulics!(soil, crop;
                             lpjmlparams::LPJmLParams = lpjmlparams)
    T = eltype(soil_management_prognostic(soil).tillage_density_factor)
    launch_1D!(
        tillage_hydraulics_kernel!,
        soil_management_prognostic(soil).tillage_density_factor,
        crop_events(crop).sowing,
        T(lpjmlparams.mixing_efficiency),
    )
    return nothing
end

function tillage_hydraulics!(state::ModelState;
                             lpjmlparams::LPJmLParams = lpjmlparams)
    density_factor = state.prognostic.soil.management.tillage_density_factor
    T = eltype(density_factor)
    launch_1D!(
        tillage_hydraulics_kernel!, density_factor, state.events.crop.sowing,
        T(lpjmlparams.mixing_efficiency),
    )
    return nothing
end

@kernel inbounds = true function tillage_hydraulics_kernel!(
    density_factor::AbstractMatrix{T},
    sowing_event::AbstractVector{S},
    mixing_efficiency::T,
) where {T <: AbstractFloat, S <: Integer}
    cell = @index(Global)
    if sowing_event[cell] != 0
        density_factor[1, cell] -=
            (density_factor[1, cell] - T(0.667)) * mixing_efficiency
    end
end
