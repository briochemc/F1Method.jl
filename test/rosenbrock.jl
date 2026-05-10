module SubModuleRosenBrock

using Test
using F1Method
using LinearAlgebra
using SciMLBase
using ForwardDiff

# Set up:
# - overload `SteadyStateProblem` solve via a tiny Newton method
# - define solver algorithm (basic Newton here)
# - define type for that algorithm (here `MyAlg`)

function newton_solve(F, ∇ₓF, x; Ftol=1e-10)
    while norm(F(x)) ≥ Ftol
        x .-= ∇ₓF(x) \ F(x)
    end
    return x
end

# Create a type for the solver's algorithm
struct MyAlg <: SciMLBase.AbstractSteadyStateAlgorithm end

# Overload SciMLBase's solve function
function SciMLBase.solve(prob::SciMLBase.AbstractSteadyStateProblem,
                          alg::MyAlg;
                          Ftol=1e-10)
    p = prob.p
    x0 = copy(prob.u0)
    F(x) = prob.f(x, p)
    ∇ₓF(x) = prob.f.jac(x, p)
    x_steady = newton_solve(F, ∇ₓF, x0, Ftol=Ftol)
    resid = F(x_steady)
    return SciMLBase.build_solution(prob, alg, x_steady, resid; retcode = ReturnCode.Success)
end


# Quasi-Rosenbrock: F(x,p) is almost the gradient of the Rosenbrock function,
# except the Hessian matrix is made non-symmetric by slightly altering F(x,p).
# With p = [a, b] the parameters defining the constants of Rosenbrock's
# function, and the slight alteration to make the Hessian non-symmetric.
# Specifically, if r(x,y) denotes Rosenbrock,
#     F₁([x, y], p) = ∂r/∂x
#     F₂([x, y], p) = 0.5 * ∂r/∂y
statefun(x, p, t = 0) = [
    -2 * (p[1] - x[1]) - 4 * p[2] * (x[2] - x[1]^2) * x[1]
    p[2] * (x[2] - x[1]^2)
]
F = ODEFunction{false}(statefun; jac = (x, p, t = 0) -> ForwardDiff.jacobian(x -> statefun(x, p), x))

# Define mismatch f(x,p) and its state-jacobian ∇ₓf(x,p)
function state_mismatch(x)
    δ(x) = x .- 1
    return 0.5δ(x)'δ(x)
end
function parameter_mismatch(p)
    δ(p) = log.(p)
    return 0.5δ(p)'δ(p)
end
f(x, p) = state_mismatch(x) + parameter_mismatch(p)
∇ₓf(x, p) = ForwardDiff.jacobian(x -> [f(x, p)], x)   # 1×n row matrix

# Initialize the starting state and parameters
x₀ = rand(2)
p₀ = rand(2)
mem = F1Method.initialize_mem(F, ∇ₓf, x₀, p₀, MyAlg())

# Define the F-1 wrappers
F1_objective(p) = F1Method.objective(f, F, mem, p, MyAlg())
F1_gradient(p)  = F1Method.gradient(f, F, ∇ₓf, mem, p, MyAlg())
F1_hessian(p)   = F1Method.hessian(f, F, ∇ₓf, mem, p, MyAlg())

# Exact solution and objective for comparison
exact_solution(p)  = [p[1], p[1]^2]
exact_objective(p) = f(exact_solution(p), p)

@testset "objective" begin
    @test exact_objective(p₀) ≈ F1_objective(p₀)
    @test exact_objective(2p₀) ≈ F1_objective(2p₀)
end
@testset "gradient" begin
    @test ForwardDiff.gradient(exact_objective, p₀)  ≈ F1_gradient(p₀)
    @test ForwardDiff.gradient(exact_objective, 2p₀) ≈ F1_gradient(2p₀)
end
@testset "Hessian" begin
    @test ForwardDiff.hessian(exact_objective, p₀)  ≈ F1_hessian(p₀)
    @test ForwardDiff.hessian(exact_objective, 2p₀) ≈ F1_hessian(2p₀)
end

end # submodule
