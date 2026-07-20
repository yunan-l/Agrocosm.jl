# using CUDA
"""
photosynthesis_C3!(PFT, photos, crop, pet, co2, temp)

Compute C3 photosynthesis rates and related diagnostic variables.
"""
function photosynthesis_C3_reference!(PFT::PftParameters,
                            crop::Crop,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T},
                            co2::AbstractArray{T};
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            photoparams::PhotoParams = photoparams,
                            comp_vcmax = false # compute vcmax internally
) where {T <: AbstractFloat}

    @unpack b = PFT
    @unpack ko25, kc25, alphac3, theta, LAMBDA_OPT = lpjmlparams
    @unpack q10ko, q10kc, po2, tau25, q10tau, cmass, cq, p, lambdamc3 = photoparams
    inactive = crop.auxiliary.photosynthesis.temperature_stress .< T(1e-2)

    ko = ko25 * q10ko .^ ((temp .- T(25.0)) * T(0.1))
    kc = kc25 * q10kc .^ ((temp .- T(25.0)) * T(0.1))
    fac = kc .* (T(1.0) .+ po2 ./ ko)
    tau = tau25 * q10tau .^ ((temp .- T(25.0)) * T(0.1)) #reflects the abiltiy of Rubisco to discriminate between CO2 and O2
    gammastar = po2 ./ (T(2.0) * tau)

    if comp_vcmax
        p_i= lambdamc3 * co2
        c1 = crop.auxiliary.photosynthesis.temperature_stress * alphac3 .* ((p_i .- gammastar) ./ (p_i .+ T(2.0) * gammastar))
        # Calculation of C2C3, Eqn 6, Haxeltine & Prentice 1996
        c2 = (p_i .- gammastar) ./ (p_i .+ fac)
        s = (24 ./ pet_daylength) * b
        sigma = 1.0f0 .- (c2 .- s) ./ (c2 .- theta * s)
        sigma = sqrt.(max.(zero(T), sigma))
        crop.auxiliary.photosynthesis.lambda .= T(LAMBDA_OPT)
        vcmax = (1.0f0 / b) * (c1 ./ c2) .* ((2.0f0 * theta - 1.0f0) .* s .- (2.0f0 * theta .* s .- c2) .* sigma) .* apar * cmass * cq
        crop.auxiliary.photosynthesis.vcmax .= ifelse.(inactive, zero(T), max.(zero(T), vcmax))
        crop.auxiliary.photosynthesis.potential_vcmax .= crop.auxiliary.photosynthesis.vcmax
        crop.auxiliary.photosynthesis.nitrogen_limitation .= ifelse.(crop.auxiliary.photosynthesis.vcmax .> zero(T), one(T), zero(T))
    end

    # calculation of C1C3, C2C3 with actual p_i (leaf internal partial pressure of CO2)
    p_i = crop.auxiliary.photosynthesis.lambda .* co2

    c1 = crop.auxiliary.photosynthesis.temperature_stress * alphac3 .* ((p_i .- gammastar) ./ (p_i .+ T(2.0) * gammastar))

    c2 = (p_i .- gammastar) ./ (p_i .+ fac)

    #   je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
    #   Convert je from daytime to hourly basis

    #   Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
    #   Eqn 3, Haxeltine & Prentice 1996

    je = c1 .* apar * cmass * cq ./ (pet_daylength .+ 1f-5)

    #   Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h
    #   Eqn 5, Haxeltine & Prentice 1996

    jc = c2 .* hour2day(crop.auxiliary.photosynthesis.vcmax)

    #   Calculation of daily gross photosynthesis, Agd, gC/m2/day
    #   Eqn 2, Haxeltine & Prentice 1996

    # round-off; a positive floor can make GPP negative at low light.
    agd = (je .+ jc .- sqrt.(max.(zero(T), (je .+ jc) .* (je .+ jc) .- T(4.0) * theta .* je .* jc))) ./ (T(2.0) * theta) .* pet_daylength
    crop.fluxes.carbon.gross_assimilation .= ifelse.(inactive, zero(T), max.(zero(T), agd))

    #   Daily dark respiration, Rd, gC/m2/day
    #   Eqn 10, Haxeltine & Prentice 1996

    #   Total daytime net photosynthesis, Adt, gC/m2/day
    #   Eqn 19, Haxeltine & Prentice 1996

    #   Daily dark respiration, Rd, gC/m2/day
    #   Eqn 10, Haxeltine & Prentice 1996
    crop.fluxes.carbon.leaf_respiration .= ifelse.(inactive, zero(T), b .* crop.auxiliary.photosynthesis.vcmax)
    adt = crop.fluxes.carbon.gross_assimilation .- hour2day(pet_daylength) .* crop.fluxes.carbon.leaf_respiration

    #   Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    crop.fluxes.carbon.net_assimilation .= max.(zero(T), adt)
    crop.fluxes.carbon.water_limited_assimilation .= ifelse.(
        adt .<= zero(T),
        zero(T),
        adt ./ cmass .* T(8.314) .* degCtoK(temp) ./ p .* T(1000.0),
    )

end


"""
photosynthesis_C4!(PFT, photos, crop, pet, co2, temp)

Compute C4 photosynthesis rates and related diagnostic variables.
"""
function photosynthesis_C4_reference!(PFT::PftParameters,
                            crop::Crop,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T};
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            photoparams::PhotoParams = photoparams,
                            comp_vcmax = false # compute vcmax internally
) where {T <: AbstractFloat}

    @unpack b = PFT
    @unpack alphac4, theta, LAMBDA_OPT = lpjmlparams
    @unpack lambdamc4, cmass, cq, p = photoparams
    inactive = crop.auxiliary.photosynthesis.temperature_stress .< T(1e-2)

    #   Parameter accounting for effect of reduced intercellular CO2
    #   concentration on photosynthesis, Phipi.
    #   Eqn 14,16, Haxeltine & Prentice 1996
    #   Fig 1b, Collatz et al 1992
    if comp_vcmax
        c1 = crop.auxiliary.photosynthesis.temperature_stress * alphac4
        c2 = 1.0f0
        s = (24 ./ pet_daylength) * b
        sigma = 1.0f0 .- (c2 .- s) ./ (c2 .- theta * s)
        # sigma = sqrt.(0.5f0 * (sigma .+ sqrt(sigma .* sigma .+ (1f-3)^2)))
        sigma = sqrt.(max.(zero(T), sigma))
        # LPJmL computes potential conductance at the common LAMBDA_OPT.
        # C4 assimilation is already saturated above lambdamc4.
        crop.auxiliary.photosynthesis.lambda .= T(LAMBDA_OPT)
        vcmax = (1.0f0 / b) * (c1 ./ c2) .* ((2.0f0 * theta - 1.0f0) .* s .- (2.0f0 * theta .* s .- c2) .* sigma) .* apar * cmass * cq
        crop.auxiliary.photosynthesis.vcmax .= ifelse.(inactive, zero(T), max.(zero(T), vcmax))
        crop.auxiliary.photosynthesis.potential_vcmax .= crop.auxiliary.photosynthesis.vcmax
        crop.auxiliary.photosynthesis.nitrogen_limitation .= ifelse.(crop.auxiliary.photosynthesis.vcmax .> zero(T), one(T), zero(T))
    end

    phipi = min.(one(T), crop.auxiliary.photosynthesis.lambda/lambdamc4)
    c1 = crop.auxiliary.photosynthesis.temperature_stress .* phipi * alphac4
    # c2 = device(ones(T, size(c1)))

    #   je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
    #   Convert je from daytime to hourly basis

    #   Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
    #   Eqn 3, Haxeltine & Prentice 1996

    je = c1 .* apar * cmass * cq ./ (pet_daylength .+ 1f-5)

    # jc = c2 .* hour2day(crop.auxiliary.photosynthesis.vcmax)
    jc = hour2day(crop.auxiliary.photosynthesis.vcmax)

    #   Calculation of daily gross photosynthesis, Agd, gC/m2/day
    #   Eqn 2, Haxeltine & Prentice 1996

    # round-off; a positive floor can make GPP negative at low light.
    agd = (je .+ jc .- sqrt.(max.(zero(T), (je .+ jc) .* (je .+ jc) .- T(4.0) * theta .* je .* jc))) ./ (T(2.0) * theta) .* pet_daylength
    crop.fluxes.carbon.gross_assimilation .= ifelse.(inactive, zero(T), max.(zero(T), agd))

    #   Daily dark respiration, Rd, gC/m2/day
    #   Eqn 10, Haxeltine & Prentice 1996

    #   Total daytime net photosynthesis, Adt, gC/m2/day
    #   Eqn 19, Haxeltine & Prentice 1996

    crop.fluxes.carbon.leaf_respiration .= ifelse.(inactive, zero(T), b .* crop.auxiliary.photosynthesis.vcmax)
    adt = crop.fluxes.carbon.gross_assimilation .- hour2day(pet_daylength) .* crop.fluxes.carbon.leaf_respiration

    #   Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    crop.fluxes.carbon.net_assimilation .= max.(zero(T), adt)
    crop.fluxes.carbon.water_limited_assimilation .= ifelse.(
        adt .<= zero(T),
        zero(T),
        adt ./ cmass .* T(8.314) .* degCtoK(temp) ./ p .* T(1000.0),
    )

end

"""Cell-local C3 photosynthesis kernel with no intermediate device arrays."""
function photosynthesis_C3!(PFT::PftParameters,
                            crop::Crop,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T},
                            co2::AbstractArray{T};
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            photoparams::PhotoParams = photoparams,
                            comp_vcmax = false
) where {T <: AbstractFloat}
    launch_1D!(
        photosynthesis_c3_kernel!,
        crop.fluxes.carbon.gross_assimilation,
        crop.fluxes.carbon.net_assimilation,
        crop.fluxes.carbon.water_limited_assimilation,
        crop.fluxes.carbon.leaf_respiration,
        crop.auxiliary.photosynthesis.potential_vcmax,
        crop.auxiliary.photosynthesis.vcmax,
        crop.auxiliary.photosynthesis.nitrogen_limitation,
        crop.auxiliary.photosynthesis.lambda,
        crop.auxiliary.photosynthesis.temperature_stress,
        apar,
        pet_daylength,
        temp,
        co2,
        PFT,
        lpjmlparams,
        photoparams,
        comp_vcmax,
    )
    return nothing
end

@kernel inbounds = true function photosynthesis_c3_kernel!(
    gross_assimilation::AbstractVector{T},
    net_assimilation::AbstractVector{T},
    water_limited_assimilation::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    potential_vcmax::AbstractVector{T},
    vcmax::AbstractVector{T},
    nitrogen_limitation::AbstractVector{T},
    lambda::AbstractVector{T},
    temperature_stress::AbstractVector{T},
    apar::AbstractVector{T},
    daylength::AbstractVector{T},
    temperature::AbstractVector{T},
    co2::AbstractVector{T},
    PFT::PftParameters,
    lpjmlparams::LPJmLParams,
    photoparams::PhotoParams,
    comp_vcmax::Bool,
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack b = PFT
    @unpack ko25, kc25, alphac3, theta, LAMBDA_OPT = lpjmlparams
    @unpack q10ko, q10kc, po2, tau25, q10tau, cmass, cq, p, lambdamc3 = photoparams

    stress = temperature_stress[cell]
    inactive = stress < T(1e-2)
    temperature_cell = temperature[cell]
    co2_cell = co2[length(co2) == 1 ? 1 : cell]
    ko = T(ko25) * T(q10ko)^((temperature_cell - T(25)) * T(0.1))
    kc = T(kc25) * T(q10kc)^((temperature_cell - T(25)) * T(0.1))
    fac = kc * (one(T) + T(po2) / ko)
    tau = T(tau25) * T(q10tau)^((temperature_cell - T(25)) * T(0.1))
    gammastar = T(po2) / (T(2) * tau)

    if comp_vcmax
        internal_co2 = T(lambdamc3) * co2_cell
        c1 = stress * T(alphac3) *
            ((internal_co2 - gammastar) / (internal_co2 + T(2) * gammastar))
        c2 = (internal_co2 - gammastar) / (internal_co2 + fac)
        s = T(24) / daylength[cell] * T(b)
        sigma = one(T) - (c2 - s) / (c2 - T(theta) * s)
        sigma = sqrt(max(zero(T), sigma))
        lambda[cell] = T(LAMBDA_OPT)
        potential = (one(T) / T(b)) * (c1 / c2) *
            ((T(2) * T(theta) - one(T)) * s -
             (T(2) * T(theta) * s - c2) * sigma) *
            apar[cell] * T(cmass) * T(cq)
        vcmax[cell] = inactive ? zero(T) : max(zero(T), potential)
        potential_vcmax[cell] = vcmax[cell]
        nitrogen_limitation[cell] = vcmax[cell] > zero(T) ? one(T) : zero(T)
    end

    internal_co2 = lambda[cell] * co2_cell
    c1 = stress * T(alphac3) *
        ((internal_co2 - gammastar) / (internal_co2 + T(2) * gammastar))
    c2 = (internal_co2 - gammastar) / (internal_co2 + fac)
    je = c1 * apar[cell] * T(cmass) * T(cq) / (daylength[cell] + T(1e-5))
    jc = c2 * hour2day(vcmax[cell])
    discriminant = max(
        zero(T),
        (je + jc) * (je + jc) - T(4) * T(theta) * je * jc,
    )
    agd = (je + jc - sqrt(discriminant)) / (T(2) * T(theta)) * daylength[cell]
    gross = inactive ? zero(T) : max(zero(T), agd)
    gross_assimilation[cell] = gross
    leaf = inactive ? zero(T) : T(b) * vcmax[cell]
    leaf_respiration[cell] = leaf
    adt = gross - hour2day(daylength[cell]) * leaf
    net_assimilation[cell] = max(zero(T), adt)
    water_limited_assimilation[cell] = adt <= zero(T) ? zero(T) :
        adt / T(cmass) * T(8.314) * (temperature_cell + T(273.15)) /
        T(p) * T(1000)
end

"""Cell-local C4 photosynthesis kernel with no intermediate device arrays."""
function photosynthesis_C4!(PFT::PftParameters,
                            crop::Crop,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T};
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            photoparams::PhotoParams = photoparams,
                            comp_vcmax = false
) where {T <: AbstractFloat}
    launch_1D!(
        photosynthesis_c4_kernel!,
        crop.fluxes.carbon.gross_assimilation,
        crop.fluxes.carbon.net_assimilation,
        crop.fluxes.carbon.water_limited_assimilation,
        crop.fluxes.carbon.leaf_respiration,
        crop.auxiliary.photosynthesis.potential_vcmax,
        crop.auxiliary.photosynthesis.vcmax,
        crop.auxiliary.photosynthesis.nitrogen_limitation,
        crop.auxiliary.photosynthesis.lambda,
        crop.auxiliary.photosynthesis.temperature_stress,
        apar,
        pet_daylength,
        temp,
        PFT,
        lpjmlparams,
        photoparams,
        comp_vcmax,
    )
    return nothing
end

@kernel inbounds = true function photosynthesis_c4_kernel!(
    gross_assimilation::AbstractVector{T},
    net_assimilation::AbstractVector{T},
    water_limited_assimilation::AbstractVector{T},
    leaf_respiration::AbstractVector{T},
    potential_vcmax::AbstractVector{T},
    vcmax::AbstractVector{T},
    nitrogen_limitation::AbstractVector{T},
    lambda::AbstractVector{T},
    temperature_stress::AbstractVector{T},
    apar::AbstractVector{T},
    daylength::AbstractVector{T},
    temperature::AbstractVector{T},
    PFT::PftParameters,
    lpjmlparams::LPJmLParams,
    photoparams::PhotoParams,
    comp_vcmax::Bool,
) where {T <: AbstractFloat}
    cell = @index(Global)
    @unpack b = PFT
    @unpack alphac4, theta, LAMBDA_OPT = lpjmlparams
    @unpack lambdamc4, cmass, cq, p = photoparams

    stress = temperature_stress[cell]
    inactive = stress < T(1e-2)
    if comp_vcmax
        c1 = stress * T(alphac4)
        s = T(24) / daylength[cell] * T(b)
        sigma = one(T) - (one(T) - s) / (one(T) - T(theta) * s)
        sigma = sqrt(max(zero(T), sigma))
        lambda[cell] = T(LAMBDA_OPT)
        potential = (one(T) / T(b)) * c1 *
            ((T(2) * T(theta) - one(T)) * s -
             (T(2) * T(theta) * s - one(T)) * sigma) *
            apar[cell] * T(cmass) * T(cq)
        vcmax[cell] = inactive ? zero(T) : max(zero(T), potential)
        potential_vcmax[cell] = vcmax[cell]
        nitrogen_limitation[cell] = vcmax[cell] > zero(T) ? one(T) : zero(T)
    end

    phipi = min(one(T), lambda[cell] / T(lambdamc4))
    c1 = stress * phipi * T(alphac4)
    je = c1 * apar[cell] * T(cmass) * T(cq) / (daylength[cell] + T(1e-5))
    jc = hour2day(vcmax[cell])
    discriminant = max(
        zero(T),
        (je + jc) * (je + jc) - T(4) * T(theta) * je * jc,
    )
    agd = (je + jc - sqrt(discriminant)) / (T(2) * T(theta)) * daylength[cell]
    gross = inactive ? zero(T) : max(zero(T), agd)
    gross_assimilation[cell] = gross
    leaf = inactive ? zero(T) : T(b) * vcmax[cell]
    leaf_respiration[cell] = leaf
    adt = gross - hour2day(daylength[cell]) * leaf
    net_assimilation[cell] = max(zero(T), adt)
    water_limited_assimilation[cell] = adt <= zero(T) ? zero(T) :
        adt / T(cmass) * T(8.314) * (temperature[cell] + T(273.15)) /
        T(p) * T(1000)
end
