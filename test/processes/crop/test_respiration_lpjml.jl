using Agrocosm
using Test

@testset "Crop NPP subtracts complete maintenance and growth respiration" begin
    for T in (Float32, Float64)
        crop = init_crop(T, 3, identity)
        crop.state.phenology.is_growing .= Int32[1, 1, 0]
        crop.state.carbon.root .= T[10, 10, 10]
        crop.state.carbon.storage .= T[4, 4, 4]
        crop.state.carbon.pool .= T[2, 2, 2]

        air = T[-20, 50, 20]
        soil = reshape(T[5, 10, 15], 1, :)
        gross = T[0, 12, 12]
        leaf = T[1, 2, 2]
        crop.fluxes.carbon.gross_assimilation .= gross
        crop.fluxes.carbon.leaf_respiration .= leaf
        respiration!(crop, cft1, air, soil, gross, leaf)

        p = lpjmlparams
        temperature_response(temp) = temp < T(-15) ? zero(T) : exp(
            T(p.e0) * (one(T) / (T(p.temp_response) + T(10)) -
                       one(T) / (min(temp, T(40)) + T(p.temp_response))),
        )
        expected_total = zeros(T, 3)
        expected_npp = zeros(T, 3)
        for cell in 1:3
            g_air = temperature_response(air[cell])
            g_soil = temperature_response(soil[1, cell])
            root = T(10) * T(cft1.respcoeff) * T(p.k) * T(cft1.nc_ratio.root) * g_soil
            storage = T(4) * T(cft1.respcoeff) * T(p.k) * T(cft1.nc_ratio.sto) * g_air
            pool = T(2) * T(cft1.respcoeff) * T(p.k) * T(cft1.nc_ratio.pool) * g_air
            growth = max(zero(T),
                         (gross[cell] - leaf[cell] - root - storage - pool) * T(p.r_growth))
            active = T(crop.state.phenology.is_growing[cell])
            expected_total[cell] = (root + storage + pool + growth) * active
            expected_npp[cell] = active == one(T) ?
                gross[cell] - leaf[cell] - expected_total[cell] : zero(T)
        end

        @test crop.fluxes.carbon.respiration ≈ expected_total rtol = T(8e-6)
        @test crop.fluxes.carbon.respiration[3] == zero(T)
        carbon_allocation!(cft1, crop)
        @test crop.fluxes.carbon.npp ≈ expected_npp rtol = T(8e-6)
    end
end
