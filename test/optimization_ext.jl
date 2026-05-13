module SubModuleOptimizationExt

using Test
using F1Method
using LinearAlgebra
using SciMLBase
using ForwardDiff
using Optimization
using OptimizationOptimJL

# Reuse the same quasi-Rosenbrock setup as test/rosenbrock.jl, but drive it
# through Optimization.jl via the F1MethodOptimizationExt extension.

function newton_solve(F, ∇ₓF, x; Ftol = 1e-10)
    while norm(F(x)) ≥ Ftol
        x .-= ∇ₓF(x) \ F(x)
    end
    return x
end

struct MyAlg <: SciMLBase.AbstractSteadyStateAlgorithm end

function SciMLBase.solve(prob::SciMLBase.AbstractSteadyStateProblem,
                          alg::MyAlg;
                          Ftol = 1e-10)
    p = prob.p
    x0 = copy(prob.u0)
    F(x) = prob.f(x, p)
    ∇ₓF(x) = prob.f.jac(x, p)
    x_steady = newton_solve(F, ∇ₓF, x0; Ftol = Ftol)
    resid = F(x_steady)
    return SciMLBase.build_solution(prob, alg, x_steady, resid; retcode = ReturnCode.Success)
end

statefun(x, p, t = 0) = [
    -2 * (p[1] - x[1]) - 4 * p[2] * (x[2] - x[1]^2) * x[1]
    p[2] * (x[2] - x[1]^2)
]
F = ODEFunction{false}(statefun; jac = (x, p, t = 0) -> ForwardDiff.jacobian(x -> statefun(x, p), x))

state_mismatch(x) = 0.5 * sum(abs2, x .- 1)
parameter_mismatch(p) = 0.5 * sum(abs2, log.(p))
f(x, p) = state_mismatch(x) + parameter_mismatch(p)
∇ₓf(x, p) = ForwardDiff.jacobian(x -> [f(x, p)], x)

x₀ = ones(2) .+ 0.1
p₀ = ones(2) .* 1.5
cache = F1Method.F1Cache(F, ∇ₓf, x₀, p₀, MyAlg())

# Build the OptimizationFunction via the extension method
optfn = F1Method.optimization_function(f, F, ∇ₓf, cache, MyAlg())

@testset "Optimization.jl extension" begin
    @test optfn isa SciMLBase.OptimizationFunction

    prob = OptimizationProblem(optfn, copy(p₀))
    res = solve(prob, NewtonTrustRegion(); maxiters = 50)

    # The objective is `state_mismatch ∘ exact_solution + parameter_mismatch`
    # = 0 + 0.5 * sum(log.(p).^2), minimised at p = [1, 1].
    @test res.u ≈ [1.0, 1.0] atol = 1e-4
end

end # submodule
