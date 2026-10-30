using NCDatasets
using Test

@testset "NetCDF output requires explicit grid mapping" begin
    mktempdir() do directory
        grid_path = joinpath(directory, "grid.nc")
        output_path = joinpath(directory, "output.nc")
        NCDataset(grid_path, "c") do ds
            defDim(ds, "longitude", 2)
            defDim(ds, "latitude", 2)
            defVar(ds, "longitude", Float32, ("longitude",))[:] = Float32[10, 11]
            defVar(ds, "latitude", Float32, ("latitude",))[:] = Float32[20, 21]
            defVar(ds, "cellid", Int32, ("longitude", "latitude"))[:, :] = Int32[1 0; 0 2]
        end

        values = Float32[1 2; 3 4]
        indices = [(1, 1), (2, 2)]
        write_output_nc(output_path, "yield", values; grid_path, grid_indices = indices)
        NCDataset(output_path) do ds
            @test ds["time"][:] == Float32[1, 2]
            @test ds["yield"][1, 1, :] == Float32[1, 3]
            @test ds["yield"][2, 2, :] == Float32[2, 4]
        end

        original = read(output_path)
        @test_throws DimensionMismatch write_output_nc(
            output_path, "yield", values;
            grid_path, grid_indices = indices, time = Float32[1],
        )
        @test_throws ArgumentError write_output_nc(
            output_path, "yield", values;
            grid_path, grid_indices = [(1, 1), (3, 2)],
        )
        @test read(output_path) == original
    end
end
