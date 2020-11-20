# For CI, make sure the downloads do not hang
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

using Test

@testset "quasi-Rosenbrock derivative" begin
    include("rosenbrock.jl")
end

@testset "AIBECS test" begin
    include("AIBECS_test.jl")
end

