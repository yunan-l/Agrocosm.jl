using Agrocosm
using Test

@testset "Pedotransfer kernel matches vector reference" begin
    cells = 7
    reference = init_soil(cells, Float32[200, 300, 500, 700, 1000], identity)
    kernel = init_soil(cells, Float32[200, 300, 500, 700, 1000], identity)

    sand = reshape(Float32.(range(0.15, 0.75; length = cells)), 1, :)
    clay = reshape(Float32.(range(0.45, 0.10; length = cells)), 1, :)
    fast = reshape(Float32.(range(40, 220; length = 5 * cells)), 5, cells)
    slow = reshape(Float32.(range(200, 900; length = 5 * cells)), 5, cells)
    saturation = reshape(Float32.(range(0.35, 0.55; length = 5 * cells)), 5, cells)
    storage = reshape(Float32.(range(20, 420; length = 5 * cells)), 5, cells)
    ice = reshape(Float32.(range(0, 40; length = 5 * cells)), 5, cells)
    density_factor = reshape(Float32.(range(0.7, 1.0; length = cells)), 1, :)
    for soil in (reference, kernel)
        soil.properties.sand_fraction .= sand
        soil.properties.clay_fraction .= clay
        soil.carbon.fast .= fast
        soil.carbon.slow .= slow
        soil.water.saturation_fraction .= saturation
        soil.water.storage .= storage
        soil.water.ice_storage .= ice
        soil.management.tillage_density_factor .= density_factor
    end

    Agrocosm.pedotransfer_reference!(reference)
    pedotransfer!(kernel)
    for field in (
        :wilting_fraction, :wilting_storage, :field_capacity,
        :saturation_fraction, :saturation_storage, :beta,
        :holding_capacity_fraction, :holding_capacity_storage,
        :saturated_conductivity, :storage, :ice_storage,
        :wilting_ice_fraction, :available_ice_storage, :free_ice_storage,
        :relative_content, :free_water,
    )
        @test isapprox(
            getproperty(kernel.water, field), getproperty(reference.water, field);
            rtol = 6.0f-6, atol = 5.0f-6,
        )
    end
end


@testset "Litter routing kernels match vector references" begin
    cells = 5
    base = init_soil(cells, soilparams.soildepth, identity)
    crop = init_crop(cells, identity)
    base.management.tillage_fraction .= Float32[
        0.05 0 0
        0.95 1 0
        0 0 1
    ]
    base.carbon.litter .= reshape(Float32.(range(1, 15; length = 3 * cells)), 3, cells)
    base.nitrogen.litter .= base.carbon.litter ./ 20.0f0
    base.carbon.input .= 0.4f0
    base.nitrogen.input .= 0.04f0

    reference = deepcopy(base)
    kernel = deepcopy(base)
    crop.events.sowing .= Int32[1, 0, 1, 0, 0]
    Agrocosm.tillage_hydraulics_reference!(reference, crop)
    tillage_hydraulics!(kernel, crop)
    @test kernel.management.tillage_density_factor ≈
        reference.management.tillage_density_factor
    Agrocosm.litter_tillage_reference!(reference, crop)
    litter_tillage!(kernel, crop)
    @test kernel.carbon.litter ≈ reference.carbon.litter
    @test kernel.nitrogen.litter ≈ reference.nitrogen.litter
    @test kernel.management.tillage_carbon ≈ reference.management.tillage_carbon
    @test kernel.management.tillage_nitrogen ≈ reference.management.tillage_nitrogen

    reference = deepcopy(base)
    kernel = deepcopy(base)
    Agrocosm.litter_bioturbation_reference!(reference)
    litter_bioturbation!(kernel)
    @test kernel.carbon.litter ≈ reference.carbon.litter
    @test kernel.nitrogen.litter ≈ reference.nitrogen.litter

    reference = deepcopy(base)
    kernel = deepcopy(base)
    crop.events.harvest .= Int32[0, 1, 0, 1, 0]
    Agrocosm.route_harvest_carbon_input_reference!(reference, crop)
    Agrocosm.route_harvest_nitrogen_input_reference!(reference, crop)
    Agrocosm.route_harvest_carbon_input!(kernel, crop)
    Agrocosm.route_harvest_nitrogen_input!(kernel, crop)
    @test kernel.carbon.litter ≈ reference.carbon.litter rtol = 3.0f-6
    @test kernel.nitrogen.litter ≈ reference.nitrogen.litter rtol = 3.0f-6
    @test kernel.management.tillage_carbon ≈ reference.management.tillage_carbon
    @test kernel.management.tillage_nitrogen ≈ reference.management.tillage_nitrogen
end

function decomposition_fixture(cells = 5)
    soil = init_soil(cells, Float32[200, 300, 500, 700, 1000], identity)
    crop = init_crop(cells, identity)
    soil.water.saturation_storage .= 180.0f0
    soil.water.holding_capacity_storage .= 100.0f0
    soil.water.wilting_storage .= 30.0f0
    soil.water.relative_content .= reshape(
        Float32.(range(0.15, 0.85; length = 5 * cells)), 5, cells,
    )
    soil.water.free_water .= 2.0f0
    soil.thermal.temperature .= reshape(
        Float32.(range(-10, 32; length = 5 * cells)), 5, cells,
    )
    soil.surface_litter.temperature .= Float32.(range(-5, 25; length = cells))
    soil.surface_litter.water_capacity .= 4.0f0
    soil.surface_litter.water_storage .= Float32.(range(0, 4; length = cells))
    soil.carbon.litter .= reshape(Float32.(range(2, 20; length = 3 * cells)), 3, cells)
    soil.nitrogen.litter .= soil.carbon.litter ./ 25.0f0
    soil.carbon.fast .= reshape(Float32.(range(20, 100; length = 5 * cells)), 5, cells)
    soil.carbon.slow .= reshape(Float32.(range(100, 500; length = 5 * cells)), 5, cells)
    soil.nitrogen.fast .= soil.carbon.fast ./ 12.0f0
    soil.nitrogen.slow .= soil.carbon.slow ./ 12.0f0
    soil.carbon.litter_response .= Float32[0.0012, 0.0008, 0.0005]
    soil.nitrogen.litter_response .= Float32[0.0012, 0.0008, 0.0005]
    vertical = Float32[0.35, 0.25, 0.18, 0.13, 0.09]
    soil.decomposition.shift_fast .= reshape(vertical, :, 1)
    soil.decomposition.shift_slow .= reshape(vertical, :, 1)
    soil.nitrogen.ammonium .= 0.4f0
    soil.nitrogen.nitrate .= 0.8f0
    soil.properties.ph .= Float32.(range(5.5, 7.5; length = cells))
    soil.management.tillage_fraction .= Float32[
        0.05 0 0
        0.95 1 0
        0 0 1
    ]
    crop.events.harvest .= Int32[0, 1, 0, 1, 0]
    soil.carbon.input .= 0.1f0
    soil.nitrogen.input .= 0.01f0
    return soil, crop
end

@testset "Soil response and C/N kernels match vector references" begin
    base, crop = decomposition_fixture()
    reference = deepcopy(base)
    kernel = deepcopy(base)
    reference_crop = deepcopy(crop)
    kernel_crop = deepcopy(crop)
    air_temperature = Float32.(range(5, 25; length = 5))
    wind = Float32.(range(1, 5; length = 5))

    Agrocosm.soil_carbon_reference!(reference_crop, reference)
    soil_carbon!(kernel_crop, kernel)
    for field in (
        :litter, :fast, :slow, :decomposed_litter, :decomposed_fast,
        :decomposed_slow, :litter_to_fast, :litter_to_slow,
        :heterotrophic_respiration,
    )
        @test isapprox(
            getproperty(kernel.carbon, field), getproperty(reference.carbon, field);
            rtol = 8.0f-6, atol = 3.0f-6,
        )
    end
    @test kernel.decomposition.response ≈ reference.decomposition.response rtol = 5.0f-6
    @test kernel.decomposition.litter_response ≈
        reference.decomposition.litter_response rtol = 5.0f-6

    Agrocosm.soil_nitrogen_reference!(
        reference_crop, reference;
        air_temperature = air_temperature, wind_speed = wind,
    )
    soil_nitrogen!(
        kernel_crop, kernel;
        air_temperature = air_temperature, wind_speed = wind,
    )
    for field in (
        :litter, :fast, :slow, :decomposed_litter, :decomposed_fast,
        :decomposed_slow, :litter_to_fast, :litter_to_slow,
        :ammonium, :nitrate, :mineralization, :immobilization,
        :nitrification, :denitrification, :volatilization,
    )
        @test isapprox(
            getproperty(kernel.nitrogen, field), getproperty(reference.nitrogen, field);
            rtol = 8.0f-6, atol = 3.0f-6,
        )
    end
end
