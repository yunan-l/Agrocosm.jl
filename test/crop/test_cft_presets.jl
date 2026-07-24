using Agrocosm
using Test

@testset "CFT presets → crop processes" begin
    @testset "phenology from the registry" begin
        wheat = CropPhenology(Float64, crop_pft(1))       # temperate cereals: laimax=7
        @test wheat.laimax ≈ 7.0
        @test wheat.fphusen ≈ 0.70 atol = 1e-6   # registry is Float32; allow the widening rounding
        maize = CropPhenology(Float64, crop_pft(3))       # maize: laimax=5
        @test maize.laimax ≈ 5.0
    end

    @testset "photosynthesis pathway + thresholds" begin
        wheat = CropPhotosynthesis(Float64, crop_pft(1))  # path 1 → C3
        @test wheat.pathway isa C3Pathway
        @test wheat.T_CO2_low ≈ 0.0 && wheat.T_CO2_high ≈ 40.0
        @test wheat.T_photos_low ≈ 12.0 && wheat.T_photos_high ≈ 17.0
        maize = CropPhotosynthesis(Float64, crop_pft(3))  # path 2 → C4
        @test maize.pathway isa C4Pathway
        @test maize.T_photos_low ≈ 21.0 && maize.T_photos_high ≈ 26.0
    end

    @testset "root distribution + nitrogen from the registry" begin
        # beta_root is CFT-specific: rice (CFT 2) differs from the others.
        @test CropRootDistribution(Float64, crop_pft(2)).beta_root ≈ 0.91 atol = 1e-6
        @test CropRootDistribution(Float64, crop_pft(1)).beta_root ≈ 0.94 atol = 1e-6
        # storage-organ C:N ratio is CFT-specific (storage_ratio): temperate roots (CFT 6) is larger.
        @test CropNitrogen(Float64, crop_pft(1)).allocation.ratio_storage ≈ 0.99 atol = 1e-6
        @test CropNitrogen(Float64, crop_pft(6)).allocation.ratio_storage ≈ 1.74 atol = 1e-6
    end

    @testset "vegetation model for a CFT" begin
        veg = CropVegetation(Float64, crop_pft("maize"))
        @test veg isa CropVegetation
        @test veg.photosynthesis.pathway isa C4Pathway
        @test veg.phenology.laimax ≈ 5.0
        @test veg.phenology_dynamics.base_temperature ≈ 5.0   # maize basetemp
        @test veg.root_distribution.beta_root ≈ 0.94 atol = 1e-6
        # a C3 crop for contrast
        wheat_veg = CropVegetation(Float64, crop_pft("temperate cereals"))
        @test wheat_veg.photosynthesis.pathway isa C3Pathway
    end
end
