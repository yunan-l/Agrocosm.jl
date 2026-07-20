function c3_adtmm_scalar_impl(lambda::T,
                              vmax::T,
                              tstress::T,
                              b::T,
                              co2::T,
                              temp::T,
                              apar::T,
                              daylength::T,
                              lpjmlparams::LPJmLParams,
                              photoparams::PhotoParams) where {T <: AbstractFloat}
    if tstress < T(1e-2)
        return zero(T)
    end

    @unpack ko25, kc25, alphac3, theta = lpjmlparams
    @unpack q10ko, q10kc, po2, tau25, q10tau, cmass, cq, p = photoparams

    ko = ko25 * q10ko^((temp - T(25)) * T(0.1))
    kc = kc25 * q10kc^((temp - T(25)) * T(0.1))
    fac = kc * (one(T) + po2 / ko)
    tau = tau25 * q10tau^((temp - T(25)) * T(0.1))
    gammastar = po2 / (T(2) * tau)
    internal_co2 = lambda * co2
    c1 = tstress * alphac3 *
         ((internal_co2 - gammastar) / (internal_co2 + T(2) * gammastar))
    c2 = (internal_co2 - gammastar) / (internal_co2 + fac)

    je = c1 * apar * cmass * cq / daylength
    jc = c2 * vmax / T(24)
    agd = (je + jc - sqrt(max(zero(T), (je + jc)^2 - T(4) * theta * je * jc))) / (T(2) * theta) * daylength
    rd = b * vmax
    adt = agd - daylength / T(24) * rd

    return adt <= zero(T) ? zero(T) :
           adt / cmass * T(8.314) * (temp + T(273.15)) / p * T(1000)
end

function c3_adtmm_scalar(lambda::T,
                         vmax::T,
                         tstress::T,
                         b::T,
                         co2::T,
                         temp::T,
                         apar::T,
                         daylength::T;
                         lpjmlparams::LPJmLParams = lpjmlparams,
                         photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    return c3_adtmm_scalar_impl(
        lambda, vmax, tstress, b, co2, temp, apar, daylength,
        lpjmlparams, photoparams,
    )
end

function c4_adtmm_scalar_impl(lambda::T,
                              vmax::T,
                              tstress::T,
                              b::T,
                              temp::T,
                              apar::T,
                              daylength::T,
                              lpjmlparams::LPJmLParams,
                              photoparams::PhotoParams) where {T <: AbstractFloat}
    if tstress < T(1e-2)
        return zero(T)
    end

    @unpack alphac4, theta = lpjmlparams
    @unpack lambdamc4, cmass, cq, p = photoparams

    phipi = min(one(T), lambda / T(lambdamc4))
    c1 = tstress * phipi * T(alphac4)
    je = c1 * apar * T(cmass) * T(cq) / daylength
    jc = vmax / T(24)
    agd = (je + jc - sqrt(max(zero(T), (je + jc)^2 - T(4) * T(theta) * je * jc))) /
          (T(2) * T(theta)) * daylength
    rd = b * vmax
    adt = agd - daylength / T(24) * rd

    return adt <= zero(T) ? zero(T) :
           adt / T(cmass) * T(8.314) * (temp + T(273.15)) / T(p) * T(1000)
end

function c4_adtmm_scalar(lambda::T,
                         vmax::T,
                         tstress::T,
                         b::T,
                         temp::T,
                         apar::T,
                         daylength::T;
                         lpjmlparams::LPJmLParams = lpjmlparams,
                         photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    return c4_adtmm_scalar_impl(
        lambda, vmax, tstress, b, temp, apar, daylength,
        lpjmlparams, photoparams,
    )
end


"""
    solve_lambda_c3_lpj(fac, vmax, tstress, b, co2, temp, apar, daylength;
                        lower=0.02, upper=0.85, tolerance=0.001,
                        max_iterations=30)

Solve LPJmL's C3 water-stress equation
`fac * (1 - lambda) - adtmm(lambda) = 0` using the compatible CPU bisection
algorithm. `co2` is atmospheric CO₂ partial pressure in Pa, while `fac` is the
water-limited conductance term constructed by the water-balance routine.

Returns `(lambda, iterations, residual)`.
"""
function solve_lambda_c3_lpj(fac::T,
                             vmax::T,
                             tstress::T,
                             b::T,
                             co2::T,
                             temp::T,
                             apar::T,
                             daylength::T;
                             lower::T = T(0.02),
                             upper::T = T(0.85),
                             tolerance::T = T(0.001),
                             max_iterations::Integer = 30,
                             lpjmlparams::LPJmLParams = lpjmlparams,
                             photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    objective(lambda) = fac * (one(T) - lambda) -
                        c3_adtmm_scalar_impl(
                            lambda, vmax, tstress, b, co2, temp, apar, daylength,
                            lpjmlparams, photoparams,
                        )

    lambda, iterations = lpj_bisect(
        objective,
        lower,
        upper;
        x_accuracy = zero(T),
        y_accuracy = tolerance,
        max_iterations = max_iterations,
    )

    return lambda, iterations, objective(lambda)
end

"""
    solve_lambda_c4_lpj(fac, vmax, tstress, b, temp, apar, daylength)

CPU reference for LPJmL's C4 water-stress lambda equation. Returns
`(lambda, iterations, residual)`.
"""
function solve_lambda_c4_lpj(fac::T,
                             vmax::T,
                             tstress::T,
                             b::T,
                             temp::T,
                             apar::T,
                             daylength::T;
                             lower::T = T(0.02),
                             upper::T = T(0.85),
                             tolerance::T = T(0.001),
                             max_iterations::Integer = 30,
                             lpjmlparams::LPJmLParams = lpjmlparams,
                             photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    objective(lambda) = fac * (one(T) - lambda) -
                        c4_adtmm_scalar_impl(
                            lambda, vmax, tstress, b, temp, apar, daylength,
                            lpjmlparams, photoparams,
                        )

    lambda, iterations = lpj_bisect(
        objective,
        lower,
        upper;
        x_accuracy = zero(T),
        y_accuracy = tolerance,
        max_iterations = max_iterations,
    )

    return lambda, iterations, objective(lambda)
end

"""
    solve_lambda_c3!(PFT, photos, crop, pet, temp, co2)

Solve the LPJmL water-stress equation independently for every grid cell. The
fixed 30-step loop and scalar, allocation-free objective are compatible with
both CPU and GPU backends. `co2` must be atmospheric partial pressure in Pa and
`crop.water.canopy_conductance` must contain actual canopy conductance after water limitation.
"""
function solve_lambda_c3!(PFT::PftParameters,
                          photos::CropPhotosynthesis,
                          crop::Crop,
                          pet::PetPar,
                          temp::AbstractArray{T},
                          co2::AbstractArray{T};
                          lpjmlparams::LPJmLParams = lpjmlparams,
                          photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    kernel_params = (
        b = T(PFT.b),
        gmin = T(PFT.gmin),
        lpjmlparams = lpjmlparams,
        photoparams = photoparams,
    )

    launch_1D!(
        solve_lambda_c3_kernel!,
        photos.lambda,
        photos.vmax,
        photos.temperature_stress,
        crop.water.canopy_conductance,
        crop.canopy.fpar,
        crop.canopy.apar,
        pet.daylength,
        temp,
        co2,
        kernel_params,
    )
end


"""
    solve_lambda_c4!(PFT, photos, crop, pet, temp, co2)

GPU/CPU backend implementation of LPJmL's C4 water-stress lambda solve.
`co2` is atmospheric partial pressure in Pa.
"""
function solve_lambda_c4!(PFT::PftParameters,
                          photos::CropPhotosynthesis,
                          crop::Crop,
                          pet::PetPar,
                          temp::AbstractArray{T},
                          co2::AbstractArray{T};
                          lpjmlparams::LPJmLParams = lpjmlparams,
                          photoparams::PhotoParams = photoparams) where {T <: AbstractFloat}
    kernel_params = (
        b = T(PFT.b),
        gmin = T(PFT.gmin),
        lpjmlparams = lpjmlparams,
        photoparams = photoparams,
    )

    launch_1D!(
        solve_lambda_c4_kernel!,
        photos.lambda,
        photos.vmax,
        photos.temperature_stress,
        crop.water.canopy_conductance,
        crop.canopy.fpar,
        crop.canopy.apar,
        pet.daylength,
        temp,
        co2,
        kernel_params,
    )
end

@kernel inbounds = true function solve_lambda_c4_kernel!(
    lambda::AbstractArray{T},
    vmax::AbstractArray{T},
    tstress::AbstractArray{T},
    conductance::AbstractArray{T},
    fpar::AbstractArray{T},
    apar::AbstractArray{T},
    daylength::AbstractArray{T},
    temp::AbstractArray{T},
    co2::AbstractArray{T},
    kernel_params,
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack b, gmin, lpjmlparams, photoparams = kernel_params
    co2_cell = co2[length(co2) == 1 ? 1 : cell]

    gpd = daylength[cell] * T(3600) *
          (conductance[cell] - gmin * fpar[cell])

    if gpd > T(1e-5) && tstress[cell] >= T(1e-2) &&
       daylength[cell] > zero(T) && co2_cell > zero(T)
        fac = gpd / T(1.6) * co2_cell * T(1e-5)
        xlow = T(0.02)
        xhigh = T(0.85)
        ylow = fac * (one(T) - xlow) -
               c4_adtmm_scalar_impl(
                   xlow, vmax[cell], tstress[cell], b, temp[cell], apar[cell],
                   daylength[cell], lpjmlparams, photoparams,
               )
        xmin = (xlow + xhigh) * T(0.5)
        ymin = typemax(T)

        for _ in 1:30
            xmid = (xlow + xhigh) * T(0.5)
            ymid = fac * (one(T) - xmid) -
                   c4_adtmm_scalar_impl(
                       xmid, vmax[cell], tstress[cell], b, temp[cell], apar[cell],
                       daylength[cell], lpjmlparams, photoparams,
                   )
            if abs(ymid) < ymin
                ymin = abs(ymid)
                xmin = xmid
            end
            if abs(ymid) < T(0.001)
                break
            elseif ylow * ymid <= zero(T)
                xhigh = xmid
            else
                xlow = xmid
                ylow = ymid
            end
        end
        lambda[cell] = xmin
    else
        # LPJmL bypasses photosynthesis and returns zero GPP here. Lambda zero
        # reproduces that result when the vector photosynthesis routine follows.
        lambda[cell] = zero(T)
    end
end

@kernel inbounds = true function solve_lambda_c3_kernel!(
    lambda::AbstractArray{T},
    vmax::AbstractArray{T},
    tstress::AbstractArray{T},
    conductance::AbstractArray{T},
    fpar::AbstractArray{T},
    apar::AbstractArray{T},
    daylength::AbstractArray{T},
    temp::AbstractArray{T},
    co2::AbstractArray{T},
    kernel_params,
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack b, gmin, lpjmlparams, photoparams = kernel_params
    co2_cell = co2[length(co2) == 1 ? 1 : cell]

    # LPJmL receives ppm and converts it to bar. Agrocosm stores partial
    # pressure in Pa, so the equivalent conversion is Pa * 1e-5.
    gpd = daylength[cell] * T(3600) *
          (conductance[cell] - gmin * fpar[cell])

    if gpd > T(1e-5) && tstress[cell] >= T(1e-2) &&
       daylength[cell] > zero(T) && co2_cell > zero(T)
        fac = gpd / T(1.6) * co2_cell * T(1e-5)
        xlow = T(0.02)
        xhigh = T(0.85)
        ylow = fac * (one(T) - xlow) -
               c3_adtmm_scalar_impl(
                   xlow, vmax[cell], tstress[cell], b, co2_cell, temp[cell],
                   apar[cell], daylength[cell], lpjmlparams, photoparams,
               )
        xmin = (xlow + xhigh) * T(0.5)
        ymin = typemax(T)

        for _ in 1:30
            xmid = (xlow + xhigh) * T(0.5)
            ymid = fac * (one(T) - xmid) -
                   c3_adtmm_scalar_impl(
                       xmid, vmax[cell], tstress[cell], b, co2_cell, temp[cell],
                       apar[cell], daylength[cell], lpjmlparams, photoparams,
                   )
            if abs(ymid) < ymin
                ymin = abs(ymid)
                xmin = xmid
            end
            if abs(ymid) < T(0.001)
                break
            elseif ylow * ymid <= zero(T)
                xhigh = xmid
            else
                xlow = xmid
                ylow = ymid
            end
        end
        lambda[cell] = xmin
    else
        # LPJmL bypasses photosynthesis and returns zero GPP here.
        lambda[cell] = zero(T)
    end
end
