using Test

gpu_tests = String[]
for (directory, _, files) in walkdir(@__DIR__)
    for file in files
        endswith(file, "_gpu.jl") || continue
        file == basename(@__FILE__) && continue
        push!(gpu_tests, joinpath(directory, file))
    end
end
sort!(gpu_tests)

@testset "Agrocosm GPU" begin
    for path in gpu_tests
        @testset "$(relpath(path, @__DIR__))" begin
            # Isolate helpers and constants declared by standalone GPU tests.
            Base.include(Module(gensym(:AgrocosmGPUTest)), path)
        end
    end
end
