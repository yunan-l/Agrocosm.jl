"""
nitrogen_transform!(soil, c_shift_fast, c_shift_slow; lpjmlparams=lpjmlparams, k_l=0.0f0)

Apply mineralization, immobilization, and nitrification transformations in soil N pools.
"""
function nitrogen_transform!(soil::Soil;
                             lpjmlparams::LPJmLParams = lpjmlparams,
                             k_l = 0.0f0 # Parton et al., 2001 equ. 2
)

    @unpack fastfrac, atmfrac, k_soil10 = lpjmlparams

    # NO3 and N2O from mineralization of litter organic matter
    # c_shift_* already includes layer-wise redistribution from equilibrium spin-up.
    F_Nmineral = sum(soil.decom_litn, dims = 1) * atmfrac .* (fastfrac * soil.c_shift_fast + (1.0f0 - fastfrac) * soil.c_shift_slow);
    soil.NH4 .+= F_Nmineral * (1 - k_l)
    soil.NO3 .+= F_Nmineral * k_l

    # NO3 and N2O from mineralization of soil organic matter
    F_Nmineral = soil.decom_fastn + soil.decom_slown
    soil.NH4 .+= F_Nmineral * (1 - k_l)
    soil.NO3 .+= F_Nmineral * k_l

    # Immobilization consumes mineral N (NH4 + NO3) and transfers it to slow soil pools.
    # decom_sum_lit* are reduced to 1D cell vectors for 1D kernel launch.
    decom_sum_litc = vec(sum(soil.decom_litc, dims = 1))
    decom_sum_litn = vec(sum(soil.decom_litn, dims = 1))
    kernel_params_immo = (lpjmlparams = lpjmlparams, cn_ratio = 15.0f0, soil_layers = 5, k_N = 5f-3)

    launch_1D!(immobilize_kernel!,
                decom_sum_litc, 
                decom_sum_litn,
                soil.NH4,
                soil.NO3,
                soil.fastn,
                soil.slown,
                soil.c_shift_fast,
                soil.c_shift_slow,
                soil.layer_depth,
                kernel_params_immo)

    # Nitrification converts NH4 to NO3 with soil moisture/temperature modifiers.
    kernel_params_nit = (lpjmlparams = lpjmlparams, soil_layers = 5, a_nit = 0.45f0, b_nit = 1.27f0, c_nit = 0.0012f0, d_nit = 2.84f0)
    
    launch_1D!(nitrify_kernel!,
                soil.ph,
                soil.NH4,
                soil.NO3,
                soil.swc,
                soil.wsats,
                soil.temp,
                kernel_params_nit)

    #  Denitrification: NO3 -> N2O + N2.
    kernel_params_denit = (lpjmlparams = lpjmlparams, soil_layers = 5)
    launch_1D!(denitrify_kernel!,
                soil.fastc,
                soil.slowc,
                soil.w,
                soil.whcs,
                soil.wpwps,
                soil.w_fw,
                soil.wsats,
                soil.temp,
                soil.NO3,
                kernel_params_denit)

    # NH3 volatilization from top-layer NH4 (no wind forcing available yet, uses default parameter).
    launch_1D!(volatilization_kernel!,
                soil.NH4,
                soil.ph,
                soil.temp,
                soil.layer_depth,
                lpjmlparams)

end


@kernel inbounds = true function immobilize_kernel!(decom_sum_litc::AbstractArray{T},
                                    decom_sum_litn::AbstractArray{T},
                                    soil_NH4::AbstractArray{M},           
                                    soil_NO3::AbstractArray{M},
                                    soil_fastn::AbstractArray{M},
                                    soil_slown::AbstractArray{M},
                                    c_shift_fast::AbstractArray{T},
                                    c_shift_slow::AbstractArray{T},
                                    soil_layer_depth::AbstractArray{T},
                                    kernel_params_immo
) where {T <: AbstractFloat, M <: AbstractFloat}
    
    cell = @index(Global)

    @unpack lpjmlparams, cn_ratio, soil_layers, k_N = kernel_params_immo
    @unpack fastfrac, atmfrac = lpjmlparams

    # Each thread updates all soil layers for one cell.
    for l in 1:soil_layers

        N_sum = soil_NH4[l, cell] + soil_NO3[l, cell]
        if(N_sum > 0) # immobilization of N 
            n_immo = fastfrac * (1 - atmfrac) * (decom_sum_litc[cell] / cn_ratio - decom_sum_litn[cell]) * c_shift_fast[l, cell] * N_sum / soil_layer_depth[l] * 1f3 / (k_N + N_sum / soil_layer_depth[l] * 1f3)
            if(n_immo > 0)
                if(n_immo > N_sum)
                    n_immo = N_sum
                end
                soil_fastn[l, cell] += n_immo
                soil_NH4[l, cell] -= n_immo * soil_NH4[l, cell] / N_sum
                soil_NO3[l, cell] -= n_immo * soil_NO3[l, cell] / N_sum
            end
        end

        # Fast/slow litter fractions are handled separately with different shift factors.
        N_sum = soil_NH4[l, cell] + soil_NO3[l, cell]
        if(N_sum > 0) # immobilization of N 
            n_immo = (1 - fastfrac) * (1 - atmfrac) * (decom_sum_litc[cell] / cn_ratio - decom_sum_litn[cell]) * c_shift_slow[l, cell] * N_sum / soil_layer_depth[l] * 1f3 / (k_N + N_sum / soil_layer_depth[l] * 1f3)
            if(n_immo > 0)
                if(n_immo > N_sum)
                    n_immo = N_sum
                end
                soil_slown[l, cell] += n_immo
                soil_NH4[l, cell] -= n_immo * soil_NH4[l, cell] / N_sum
                soil_NO3[l, cell] -= n_immo * soil_NO3[l, cell] / N_sum
            end
        end
    end

end


@kernel inbounds = true function nitrify_kernel!(
                                 soil_ph::AbstractArray{T},
                                 soil_NH4::AbstractArray{M},           
                                 soil_NO3::AbstractArray{M},
                                 soil_swc::AbstractArray{M},
                                 soil_wsats::AbstractArray{M},
                                 soil_temp::AbstractArray{M},
                                 kernel_params_nit
) where {T <: AbstractFloat, M <: AbstractFloat}
    
    cell = @index(Global)

    @unpack lpjmlparams, soil_layers, a_nit, b_nit, c_nit, d_nit = kernel_params_nit
    @unpack k_max, k_2 = lpjmlparams

    # Potential nitrification rate is shaped by water-filled pore space and temperature response.
    for l in 1:soil_layers

        x = soil_swc[l, cell] / max(soil_wsats[l, cell], T(1e-8))
        n_nit = a_nit - b_nit
        m_nit = a_nit - c_nit
        z_nit = d_nit * (b_nit - a_nit) / (a_nit - c_nit)
        base1 = (x - b_nit) / n_nit
        base2 = (x - c_nit) / m_nit
        if base1 <= zero(T) || base2 <= zero(T)
            fac_wfps = zero(T)
        else
            fac_wfps = base1^(z_nit) * base2^(d_nit)
            if !isfinite(fac_wfps) || fac_wfps < zero(T)
                fac_wfps = zero(T)
            end
        end
        fac_temp = exp(-(soil_temp[l, cell] - T(18.79))^2 / T(2*5.26*5.26))
        fac_ph = T(0.56) + atan(T(π) * T(0.45) * (soil_ph[cell] - T(5.0))) / T(π)

        F_NO3 = k_max * soil_NH4[l, cell] * fac_temp * fac_wfps * fac_ph
        if F_NO3 > soil_NH4[l, cell]
            F_NO3 = soil_NH4[l, cell]
        end
        # F_N2O = k_2 * F_NO3
        soil_NO3[l, cell] += F_NO3 * (1 - k_2)
        soil_NH4[l, cell] -= F_NO3
    end
end


@kernel inbounds = true function denitrify_kernel!(
                                 soil_fastc::AbstractArray{M},
                                 soil_slowc::AbstractArray{M},
                                 soil_w::AbstractArray{M},
                                 soil_whcs::AbstractArray{M},
                                 soil_wpwps::AbstractArray{M},
                                 soil_w_fw::AbstractArray{M},
                                 soil_wsats::AbstractArray{M},
                                 soil_temp::AbstractArray{M},
                                 soil_NO3::AbstractArray{M},
                                 kernel_params_denit
) where {M <: AbstractFloat}

    cell = @index(Global)

    @unpack lpjmlparams, soil_layers = kernel_params_denit
    @unpack CDN, n2o_denit_frac = lpjmlparams

    for l in 1:soil_layers
        Corg = max(zero(M), soil_fastc[l, cell] + soil_slowc[l, cell])
        temp = soil_temp[l, cell]

        FT = if temp > zero(M)
            M(0.0326) + M(0.00351) * temp^M(1.652) - (temp / M(41.748))^M(7.19)
        elseif temp > M(45.9)
            zero(M)
        else
            M(0.0326)
        end
        denit_t = (soil_wpwps[l, cell] + soil_w[l, cell] * soil_whcs[l, cell] + soil_w_fw[l, cell]) / max(soil_wsats[l, cell], M(1e-8))

        N_denit = zero(M)
        if temp <= M(45.9)
            FW = min(one(M), M(6.664096e-10) * exp(M(21.12912) * denit_t))
            TCDF = one(M) - exp(-CDN * FT * Corg)
            N_denit = FW * TCDF * soil_NO3[l, cell]
        end
        N_denit = min(max(N_denit, zero(M)), soil_NO3[l, cell])
        soil_NO3[l, cell] -= N_denit

        # Keep this split for parity with LPJmL even if not emitted to outputs here.
        N2O_denit = n2o_denit_frac * N_denit
        _N2_denit = N_denit - N2O_denit
    end
end


@kernel inbounds = true function volatilization_kernel!(
                                     soil_NH4::AbstractArray{M},
                                     soil_ph::AbstractArray{T},
                                     soil_temp::AbstractArray{M},
                                     soil_layer_depth::AbstractArray{T},
                                     lpjmlparams::LPJmLParams
) where {T <: AbstractFloat, M <: AbstractFloat}

    cell = @index(Global)
    
    @unpack volatil_wind, volatil_length = lpjmlparams

    temp = soil_temp[1, cell]
    pH = soil_ph[cell]
    NH4 = max(zero(M), soil_NH4[1, cell])

    # LPJmL volatilization.c (Montes 2009 parameterization)
    k_a = M(10)^(M(0.05) - M(2788.0) / (temp + M(273.15)))
    f_nh3 = one(M) / (one(M) + M(10)^(-pH) / max(k_a, M(1e-12)))
    nh3_solution = f_nh3 * NH4 / max(soil_layer_depth[1], T(1e-8)) * M(1000.0)
    k_h = M(0.2138) / (temp + M(273.15)) * M(10)^(M(6.123) - M(1825.0) / (temp + M(273.15)))
    nh3_gas = k_h * nh3_solution
    h_m = M(0.000612) * volatil_wind^M(0.8) * (temp + M(273.15))^M(0.382) * volatil_length^M(-0.2)
    vol_flux = M(86400.0) * h_m * nh3_gas
    vol_flux = min(max(vol_flux, zero(M)), soil_NH4[1, cell])
    soil_NH4[1, cell] -= vol_flux
end
