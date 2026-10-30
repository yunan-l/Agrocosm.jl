using CUDA

CUDA.functional() || error("A functional NVIDIA GPU is required for this test")
CUDA.allowscalar(false)

include("test_harvest_balance.jl")

@testset "CUDA last-harvest annual yield" begin
    test_last_harvest_yield(CuArray)
end

@testset "CUDA failed crop returns unallocated nitrogen" begin
    cells = 1024
    crop = init_crop(cells, CuArray)
    soil = init_soil(cells, soilparams.soildepth, CuArray)
    output = init_output(cells, CuArray)
    state = test_model_state(crop, soil; output)
    Agrocosm.prepare_output_block!(output, 1, 0)
    crop.events.harvest .= Int32(1)
    crop.state.phenology.is_growing .= Int32(1)
    crop.state.nitrogen.total .= 0.5f0
    crop.state.nitrogen.leaf .= 0.1f0
    crop.state.nitrogen.root .= 0.1f0
    crop.state.nitrogen.pool .= 0.2f0

    Agrocosm.terminate_failed_crop!(
        state, state, output, CUDA.fill(0.6f0, cells), 100; output_row = 1,
    )
    synchronize()
    @test all(≈(0.2f0), Array(soil.nitrogen.input)[Agrocosm.ROOT_LITTER, :])
    @test all(≈(0.5f0), Array(crop.fluxes.nitrogen.harvest_export) .+
                         vec(sum(Array(soil.nitrogen.input); dims = 1)))
end
