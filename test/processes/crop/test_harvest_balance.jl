using Agrocosm
using Test

@testset "Harvest carbon and nitrogen routing is conservative" begin
    crop = init_crop(1, identity)
    soil = init_soil(1, soilparams.soildepth, identity)
    output = init_output(1, identity)
    residue_fraction = Float32[0.6]
    soil.management.tillage_fraction .= Float32[1 0 0; 0 1 0; 0 0 1]

    crop.state.carbon.leaf .= 2.0f0
    crop.state.carbon.root .= 4.0f0
    crop.state.carbon.storage .= 5.0f0
    crop.state.carbon.pool .= 3.0f0
    crop.state.carbon.biomass .= 14.0f0
    crop.state.nitrogen.leaf .= 0.2f0
    crop.state.nitrogen.root .= 0.4f0
    crop.state.nitrogen.storage .= 0.5f0
    crop.state.nitrogen.pool .= 0.3f0
    crop.state.nitrogen.total .= 1.4f0
    crop.state.phenology.is_growing .= Int32(1)
    crop.state.phenology.harvesting_previous .= false
    crop.state.phenology.harvesting .= true

    carbon_before = crop.state.carbon.leaf[1] + crop.state.carbon.root[1] +
        crop.state.carbon.storage[1] + crop.state.carbon.pool[1]
    nitrogen_before = crop.state.nitrogen.total[1]
    harvest_crop!(crop, soil, output, residue_fraction, 150)

    carbon_export = crop.fluxes.carbon.harvest_export[1]
    nitrogen_export = crop.fluxes.nitrogen.harvest_export[1]
    carbon_residue = sum(soil.carbon.input)
    nitrogen_residue = sum(soil.nitrogen.input)
    @test crop.events.harvest[1] == 1
    @test crop.fluxes.carbon.yield[1] == 5.0f0
    @test carbon_export == 7.0f0
    @test carbon_residue == 7.0f0
    @test carbon_export + carbon_residue == carbon_before
    @test nitrogen_export == 0.7f0
    @test nitrogen_residue ≈ 0.7f0 atol = 1.0f-7
    @test nitrogen_export + nitrogen_residue ≈ nitrogen_before atol = 2.0f-7

    Agrocosm.route_harvest_residues!(soil, crop)
    @test sum(soil.carbon.litter) == carbon_residue
    @test sum(soil.nitrogen.litter) == nitrogen_residue
    @test crop.state.phenology.is_growing[1] == 0
end
