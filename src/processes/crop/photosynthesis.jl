# using CUDA
"""
photosynthesis_C3!(PFT, photos, crop, pet, co2, temp)

Compute C3 photosynthesis rates and related diagnostic variables.
"""
function photosynthesis_C3!(PFT::PftParameters,
                            photos::Photos,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T},
                            co2::AbstractArray{T};
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            photoparams::PhotoParams = photoparams,
                            comp_vmax = false # compute vmax internally
) where {T <: AbstractFloat}
    
    @unpack b = PFT
    @unpack ko25, kc25, alphac3, theta = lpjmlparams
    @unpack q10ko, q10kc, po2, tau25, q10tau, cmass, cq, p, lambdamc3 = photoparams
    inactive = photos.tstress .< T(1e-2)
    
    ko = ko25 * q10ko .^ ((temp .- T(25.0)) * T(0.1))
    kc = kc25 * q10kc .^ ((temp .- T(25.0)) * T(0.1))
    fac = kc .* (T(1.0) .+ po2 ./ ko)
    tau = tau25 * q10tau .^ ((temp .- T(25.0)) * T(0.1)) #reflects the abiltiy of Rubisco to discriminate between CO2 and O2
    gammastar = po2 ./ (T(2.0) * tau)

    if comp_vmax
        p_i= lambdamc3 * co2
        c1 = photos.tstress * alphac3 .* ((p_i .- gammastar) ./ (p_i .+ T(2.0) * gammastar))
        # Calculation of C2C3, Eqn 6, Haxeltine & Prentice 1996
        c2 = (p_i .- gammastar) ./ (p_i .+ fac)
        s = (24 ./ pet_daylength) * b
        sigma = 1.0f0 .- (c2 .- s) ./ (c2 .- theta * s)
        sigma = sqrt.(max.(zero(T), sigma))
        photos.lambda .= 0.8f0  
        vmax = (1.0f0 / b) * (c1 ./ c2) .* ((2.0f0 * theta - 1.0f0) .* s .- (2.0f0 * theta .* s .- c2) .* sigma) .* apar * cmass * cq
        photos.vmax = ifelse.(inactive, zero(T), max.(zero(T), vmax))
    end

    # calculation of C1C3, C2C3 with actual p_i (leaf internal partial pressure of CO2)
    p_i = photos.lambda .* co2

    c1 = photos.tstress * alphac3 .* ((p_i .- gammastar) ./ (p_i .+ T(2.0) * gammastar))

    c2 = (p_i .- gammastar) ./ (p_i .+ fac)

    #   je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
    #   Convert je from daytime to hourly basis

    #   Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
    #   Eqn 3, Haxeltine & Prentice 1996

    je = c1 .* apar * cmass * cq ./ (pet_daylength .+ 1f-5)

    #   Calculation of rubisco-activity-limited photosynthesis rate JC, molC/m2/h
    #   Eqn 5, Haxeltine & Prentice 1996

    jc = c2 .* hour2day(photos.vmax)

    #   Calculation of daily gross photosynthesis, Agd, gC/m2/day
    #   Eqn 2, Haxeltine & Prentice 1996

    # round-off; a positive floor can make GPP negative at low light.
    agd = (je .+ jc .- sqrt.(max.(zero(T), (je .+ jc) .* (je .+ jc) .- T(4.0) * theta .* je .* jc))) ./ (T(2.0) * theta) .* pet_daylength
    photos.agd = ifelse.(inactive, zero(T), max.(zero(T), agd))

    #   Daily dark respiration, Rd, gC/m2/day
    #   Eqn 10, Haxeltine & Prentice 1996

    #   Total daytime net photosynthesis, Adt, gC/m2/day
    #   Eqn 19, Haxeltine & Prentice 1996

    #   Daily dark respiration, Rd, gC/m2/day
    #   Eqn 10, Haxeltine & Prentice 1996
    photos.rd .= ifelse.(inactive, zero(T), b .* photos.vmax)
    adt = photos.agd .- hour2day(pet_daylength) .* photos.rd

    #   Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    photos.adt = max.(zero(T), adt)
    photos.adtmm = ifelse.(
        adt .<= zero(T),
        zero(T),
        adt ./ cmass .* T(8.314) .* degCtoK(temp) ./ p .* T(1000.0),
    )

end


"""
photosynthesis_C4!(PFT, photos, crop, pet, co2, temp)

Compute C4 photosynthesis rates and related diagnostic variables.
"""
function photosynthesis_C4!(PFT::PftParameters,
                            photos::Photos,
                            apar::AbstractArray{T},
                            pet_daylength::AbstractArray{T},
                            temp::AbstractArray{T};
                            lpjmlparams::LPJmLParams = lpjmlparams,
                            photoparams::PhotoParams = photoparams,
                            comp_vmax = false # compute vmax internally
) where {T <: AbstractFloat}
    
    @unpack b = PFT
    @unpack alphac4, theta = lpjmlparams
    @unpack lambdamc4, cmass, cq, p = photoparams
    inactive = photos.tstress .< T(1e-2)
    
    #   Parameter accounting for effect of reduced intercellular CO2
    #   concentration on photosynthesis, Phipi.
    #   Eqn 14,16, Haxeltine & Prentice 1996
    #   Fig 1b, Collatz et al 1992
    if comp_vmax
        c1 = photos.tstress * alphac4
        c2 = 1.0f0
        s = (24 ./ pet_daylength) * b
        sigma = 1.0f0 .- (c2 .- s) ./ (c2 .- theta * s)
        # sigma = sqrt.(0.5f0 * (sigma .+ sqrt(sigma .* sigma .+ (1f-3)^2)))
        sigma = sqrt.(max.(zero(T), sigma))
        photos.lambda .= 0.4f0  
        vmax = (1.0f0 / b) * (c1 ./ c2) .* ((2.0f0 * theta - 1.0f0) .* s .- (2.0f0 * theta .* s .- c2) .* sigma) .* apar * cmass * cq
        photos.vmax = ifelse.(inactive, zero(T), max.(zero(T), vmax))
    end

    phipi = min.(one(T), photos.lambda/lambdamc4)
    c1 = photos.tstress .* phipi * alphac4
    # c2 = device(ones(T, size(c1)))

    #   je is PAR-limited photosynthesis rate molC/m2/h, Eqn 3
    #   Convert je from daytime to hourly basis

    #   Calculation of PAR-limited photosynthesis rate, JE, molC/m2/h
    #   Eqn 3, Haxeltine & Prentice 1996

    je = c1 .* apar * cmass * cq ./ (pet_daylength .+ 1f-5)
    
    # jc = c2 .* hour2day(photos.vmax)
    jc = hour2day(photos.vmax)

    #   Calculation of daily gross photosynthesis, Agd, gC/m2/day
    #   Eqn 2, Haxeltine & Prentice 1996

    # round-off; a positive floor can make GPP negative at low light.
    agd = (je .+ jc .- sqrt.(max.(zero(T), (je .+ jc) .* (je .+ jc) .- T(4.0) * theta .* je .* jc))) ./ (T(2.0) * theta) .* pet_daylength
    photos.agd = ifelse.(inactive, zero(T), max.(zero(T), agd))

    #   Daily dark respiration, Rd, gC/m2/day
    #   Eqn 10, Haxeltine & Prentice 1996

    #   Total daytime net photosynthesis, Adt, gC/m2/day
    #   Eqn 19, Haxeltine & Prentice 1996

    photos.rd .= ifelse.(inactive, zero(T), b .* photos.vmax)
    adt = photos.agd .- hour2day(pet_daylength) .* photos.rd

    #   Convert adt from gC/m2/day to mm/m2/day using ideal gas equation
    photos.adt = max.(zero(T), adt)
    photos.adtmm = ifelse.(
        adt .<= zero(T),
        zero(T),
        adt ./ cmass .* T(8.314) .* degCtoK(temp) ./ p .* T(1000.0),
    )

end
