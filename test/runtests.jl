# For CI, make sure the downloads do not hang
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

using Test

@testset "dense (quasi-Rosenbrock)" begin
    include("dense_rosenbrock.jl")
end

@testset "Optimization.jl extension" begin
    include("optimization_ext.jl")
end

@testset "sparse (AIBECS P-model)" begin
    include("sparse_aibecs.jl")
end
