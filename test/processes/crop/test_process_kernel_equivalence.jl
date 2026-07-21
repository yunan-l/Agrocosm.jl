using Agrocosm
using Test

function compare_fields(reference, kernel, fields; rtol = 2.0f-6, atol = 2.0f-7)
    for field in fields
        @test getproperty(kernel, field) ≈ getproperty(reference, field) rtol = rtol atol = atol
    end
end


@testset "Canopy radiation kernels match vector references" begin
    cells = 8
    for pft in (cft1, cft3)
        reference = init_crop(cells, identity)
        kernel = init_crop(cells, identity)
        soil_reference = init_soil(cells, soilparams.soildepth, identity)
        soil_kernel = init_soil(cells, soilparams.soildepth, identity)
        pet_reference = init_pet(cells, identity)
        pet_kernel = init_pet(cells, identity)
        phenology_fraction = Float32.(range(0, 1; length = cells))
        lai = Float32(pft.laimax) .* phenology_fraction
        growing = Int32[0, 1, 1, 0, 1, 1, 0, 1]
        par = Float32.(range(0, 25; length = cells))
        litter_carbon = Float32.(range(0, 250; length = cells))
        snow_height = Float32[0, 0, 0.1, 0.2, 0, 0.05, 0, 0.1]
        snow_fraction = Float32[0, 0, 0.4, 0.8, 0, 0.2, 0, 1]
        for (crop, pet) in ((reference, pet_reference), (kernel, pet_kernel))
            crop.state.canopy.lai .= lai
            crop.state.phenology.is_growing .= growing
            pet.par .= par
        end
        for soil in (soil_reference, soil_kernel)
            soil.carbon.litter[1, :] .= litter_carbon
            soil.snow.height .= snow_height
            soil.snow.fraction .= snow_fraction
        end
        maize = pft === cft3
        Agrocosm.albedo_reference!(
            pft, reference, soil_reference, pet_reference; maize = maize,
        )
        albedo!(pft, kernel, soil_kernel, pet_kernel; maize = maize)
        @test kernel.auxiliary.canopy.albedo ≈ reference.auxiliary.canopy.albedo rtol = 2.0f-6
        @test pet_kernel.albedo ≈ pet_reference.albedo rtol = 2.0f-6

        if pft === cft3
            Agrocosm.apar_crop_maize_reference!(pft, reference, pet_reference)
            apar_crop_maize!(pft, kernel, pet_kernel)
        else
            Agrocosm.apar_crop_reference!(pft, reference, pet_reference)
            apar_crop!(pft, kernel, pet_kernel)
        end
        @test kernel.auxiliary.canopy.fpar ≈ reference.auxiliary.canopy.fpar rtol = 2.0f-6
        @test kernel.auxiliary.canopy.apar ≈ reference.auxiliary.canopy.apar rtol = 2.0f-6
    end
end

@testset "Cultivation kernel matches vector reference" begin
    cells = 6
    reference = init_crop(cells, identity)
    kernel = init_crop(cells, identity)
    reference_soil = init_soil(cells, soilparams.soildepth, identity)
    kernel_soil = init_soil(cells, soilparams.soildepth, identity)
    reference_land = init_managed_land(cells, identity)
    kernel_land = init_managed_land(cells, identity)
    sowing_dates = Int32[100, 99, 100, 101, 100, 200]
    reference.auxiliary.calendar.sowing_date .= sowing_dates
    kernel.auxiliary.calendar.sowing_date .= sowing_dates
    reference.state.carbon.biomass .= 3.0f0
    kernel.state.carbon.biomass .= 3.0f0
    reference.state.nitrogen.total .= 0.2f0
    kernel.state.nitrogen.total .= 0.2f0
    Agrocosm.cultivate_reference!(
        reference, reference_land, reference_soil, 100;
        apply_prescribed_fertilizer = false,
    )
    cultivate!(
        kernel, kernel_land, kernel_soil, 100;
        apply_prescribed_fertilizer = false,
    )
    for (reference_container, kernel_container, fields) in (
        (reference.state.phenology, kernel.state.phenology, (:harvesting, :is_growing)),
        (reference.state.canopy, kernel.state.canopy, (:lai,)),
        (reference.state.carbon, kernel.state.carbon, (:biomass, :root, :leaf, :storage, :pool)),
        (reference.state.nitrogen, kernel.state.nitrogen, (:total,)),
        (reference.fluxes.nitrogen, kernel.fluxes.nitrogen, (:seed_input,)),
        (reference.events, kernel.events, (:sowing,)),
    )
        for field in fields
            @test getproperty(kernel_container, field) ≈ getproperty(reference_container, field)
        end
    end
end

@testset "Respiration kernel matches vector reference" begin
    cells = 8
    reference = init_crop(cells, identity)
    kernel = init_crop(cells, identity)
    root = Float32.(range(1, 30; length = cells))
    storage = Float32.(range(0, 20; length = cells))
    pool = Float32.(range(2, 12; length = cells))
    growing = Int32[1, 1, 0, 1, 0, 1, 1, 1]
    temperature = Float32[-45, -20, 0, 10, 20, 30, 35, 40]
    soil_temperature = reshape(Float32[-20, -10, 0, 5, 10, 15, 20, 25], 1, :)
    gross = Float32.(range(0, 18; length = cells))
    leaf_respiration = Float32.(range(0, 2; length = cells))
    for crop in (reference, kernel)
        crop.state.carbon.root .= root
        crop.state.carbon.storage .= storage
        crop.state.carbon.pool .= pool
        crop.state.phenology.is_growing .= growing
    end
    destination = kernel.fluxes.carbon.respiration
    Agrocosm.respiration_reference!(
        reference, cft1, temperature, soil_temperature, gross .- leaf_respiration,
    )
    respiration!(kernel, cft1, temperature, soil_temperature, gross, leaf_respiration)
    @test kernel.fluxes.carbon.respiration === destination
    @test kernel.fluxes.carbon.respiration ≈ reference.fluxes.carbon.respiration rtol = 2.0f-6 atol = 2.0f-7
end

@testset "PET/PAR kernel matches vector reference" begin
    cells = 8
    for T in (Float32, Float64)
        reference = init_pet(T, cells, identity)
        kernel = init_pet(T, cells, identity)
        albedo = T.(range(0.1, 0.4; length = cells))
        reference.albedo .= albedo
        kernel.albedo .= albedo
        latitude = T[-70, -45, -10, 0, 10, 45, 70, 80]
        temperature = T[-20, -5, 0, 10, 20, 30, 35, 40]
        longwave = T.(range(-120, 40; length = cells))
        shortwave = T.(range(0, 350; length = cells))
        daylength_destination = kernel.daylength
        Agrocosm.petpar_reference!(
            reference, 172, latitude, temperature, longwave, shortwave,
        )
        petpar!(kernel, 172, latitude, temperature, longwave, shortwave)
        @test kernel.daylength === daylength_destination
        compare_fields(reference, kernel, (:daylength, :par, :eeq); rtol = T(3e-6))
    end
end

@testset "C3/C4 kernels match vector references" begin
    cells = 8
    apar = Float32[0, 2, 5, 10, 15, 20, 25, 30]
    daylength = Float32[6, 8, 10, 12, 14, 16, 18, 20]
    temperature = Float32[-10, 0, 10, 20, 25, 30, 35, 40]
    co2 = Float32[40]
    stress = Float32[0, 0.005, 0.2, 0.5, 0.8, 1, 0.7, 0.3]
    for (pft, reference_function, kernel_function) in (
        (cft1, Agrocosm.photosynthesis_C3_reference!, photosynthesis_C3!),
        (cft3, Agrocosm.photosynthesis_C4_reference!, photosynthesis_C4!),
    )
        reference = init_crop(cells, identity)
        kernel = init_crop(cells, identity)
        reference.auxiliary.photosynthesis.temperature_stress .= stress
        kernel.auxiliary.photosynthesis.temperature_stress .= stress
        gross_destination = kernel.fluxes.carbon.gross_assimilation
        if pft === cft1
            reference_function(pft, reference, apar, daylength, temperature, co2; comp_vcmax = true)
            kernel_function(pft, kernel, apar, daylength, temperature, co2; comp_vcmax = true)
        else
            reference_function(pft, reference, apar, daylength, temperature; comp_vcmax = true)
            kernel_function(pft, kernel, apar, daylength, temperature; comp_vcmax = true)
        end
        @test kernel.fluxes.carbon.gross_assimilation === gross_destination
        compare_fields(
            reference.fluxes.carbon, kernel.fluxes.carbon,
            (:gross_assimilation, :net_assimilation, :water_limited_assimilation,
             :leaf_respiration);
            rtol = 4.0f-6, atol = 3.0f-7,
        )
        compare_fields(
            reference.auxiliary.photosynthesis, kernel.auxiliary.photosynthesis,
            (:potential_vcmax, :vcmax, :nitrogen_limitation, :lambda);
            rtol = 4.0f-6, atol = 3.0f-7,
        )

        new_lambda = Float32.(range(0.45, 0.85; length = cells))
        reference.auxiliary.photosynthesis.lambda .= new_lambda
        kernel.auxiliary.photosynthesis.lambda .= new_lambda
        if pft === cft1
            reference_function(pft, reference, apar, daylength, temperature, co2; comp_vcmax = false)
            kernel_function(pft, kernel, apar, daylength, temperature, co2; comp_vcmax = false)
        else
            reference_function(pft, reference, apar, daylength, temperature; comp_vcmax = false)
            kernel_function(pft, kernel, apar, daylength, temperature; comp_vcmax = false)
        end
        compare_fields(
            reference.fluxes.carbon, kernel.fluxes.carbon,
            (:gross_assimilation, :net_assimilation, :water_limited_assimilation,
             :leaf_respiration);
            rtol = 4.0f-6, atol = 3.0f-7,
        )
        @test reference.auxiliary.photosynthesis.vcmax ≈ kernel.auxiliary.photosynthesis.vcmax
    end
end
