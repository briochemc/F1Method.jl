# This AIBECS test is derived from the P-model tutorial in AIBECS.jl
# but uses one of the small circulations to test derivatives of the cost
# function against finite differences.
module SubModuleAIBECS

using Test
using LinearAlgebra
using F1Method
using AIBECS
import AIBECS: @units, units
import AIBECS: @limits, limits
import AIBECS: @initial_value, initial_value
import AIBECS: @flattenable, flattenable, flatten
import AIBECS: @prior, prior
using Unitful: m, d, s, yr, Myr, mol, mmol, μmol, μM
using Distributions
using Bijectors
using DataFrames
using WorldOceanAtlasTools
using FiniteDiff
using ForwardDiff
using NonlinearSolve
using LinearSolve

# AIBECS model — tiny 2x2x2 grid for fast tests
grd, T_Circ = Primeau_2x2x2.load()
T_DIP(p) = T_Circ
T_POP(p) = transportoperator(grd, z -> w(z, p))
function w(z, p)
    @unpack w₀, w′ = p
    w₀ + w′ * z
end
z_top = topdepthvec(grd) # uptake only in top layer
function U(x, p)
    @unpack σ, τ_DIP, k = p
    return @. σ * x / τ_DIP * x / (x + k) * (z_top == 0) * (x ≥ 0)
end
function R(x, p)
    @unpack τ_POP = p
    return x / τ_POP
end
function G_DIP(DIP, POP, p)
    @unpack DIP_geo, τ_geo = p
    return @. -$U(DIP, p) + $R(POP, p) + (DIP_geo - DIP) / τ_geo
end
function G_POP(DIP, POP, p)
    @unpack τ_geo = p
    return @. $U(DIP, p) - $R(POP, p) - POP / τ_geo
end

# Analytical ∂Gs (diagonal vectors of the partial derivatives of G_i wrt x_j).
# Used via the new AIBECSFunction(...; ∂Gs = ...) keyword to skip ForwardDiff
# inside the AIBECS Jacobian assembly — that lets us put ForwardDiff *around*
# the steady-state solve as a clean ground truth for F1Method.gradient and
# F1Method.hessian, with no nested-Dual contamination.
#
# d/dx[x²/(x+k)] = x*(x+2k)/(x+k)²
function ∂U_∂DIP(DIP, p)
    @unpack σ, τ_DIP, k = p
    in_euphotic = z_top .== 0
    return @. σ / τ_DIP * in_euphotic * (DIP ≥ 0) * DIP * (DIP + 2k) / (DIP + k)^2
end
function ∂G_DIP_∂DIP(DIP, POP, p)
    @unpack τ_geo = p
    return @. -$∂U_∂DIP(DIP, p) - 1 / τ_geo
end
function ∂G_DIP_∂POP(DIP, POP, p)
    @unpack τ_POP = p
    return @. 1 / τ_POP * one(POP)
end
function ∂G_POP_∂DIP(DIP, POP, p)
    return ∂U_∂DIP(DIP, p)
end
function ∂G_POP_∂POP(DIP, POP, p)
    @unpack τ_POP, τ_geo = p
    return @. -1 / τ_POP - 1 / τ_geo + 0 * POP   # 0*POP gives the right shape
end
∂Gs = (
    (∂G_DIP_∂DIP, ∂G_DIP_∂POP),
    (∂G_POP_∂DIP, ∂G_POP_∂POP),
)

∞ = Inf
@initial_value @units @flattenable @limits struct PmodelParameters{U} <: AbstractParameters{U}
    w₀::U       |  0.64 | m/d      | true  | (0, ∞)
    w′::U       |  0.13 | m/d/m    | true  | (0, ∞)
    τ_DIP::U    | 230.0 | d        | true  | (0, ∞)
    k::U        |  6.62 | μmol/m^3 | true  | (0, ∞)
    τ_POP::U    |   5.0 | d        | true  | (0, ∞)
    τ_geo::U    |   1.0 | Myr      | false | (0, ∞)
    DIP_geo::U  |  2.12 | mmol/m^3 | true  | (-∞, ∞)
    σ::U        |   0.3 | NoUnits  | true  | (0, 1)
end
function prior(::Type{T}, s::Symbol) where {T <: AbstractParameters}
    if flattenable(T, s)
        lb, ub = limits(T, s)
        if (lb, ub) == (0, ∞)
            μ = log(initial_value(T, s))
            LogNormal(μ, 1.0)
        elseif (lb, ub) == (-∞, ∞)
            μ = initial_value(T, s)
            σ = 10.0
            Distributions.Normal(μ, σ)
        else
            m = initial_value(T, s)
            f = (m - lb) / (ub - lb)
            LocationScale(lb, ub - lb, LogitNormal(log(f / (1 - f)), 1.0))
        end
    else
        nothing
    end
end
prior(::T, s::Symbol) where {T <: AbstractParameters} = prior(T, s)
prior(::Type{T}) where {T <: AbstractParameters} = Tuple(prior(T, s) for s in AIBECS.symbols(T))
prior(::T) where {T <: AbstractParameters} = prior(T)

p = PmodelParameters()
λ = p2λ(p)
nb = sum(iswet(grd))
# Use analytical ∂Gs so the inner AIBECS Jacobian is a closed-form sparse
# matrix instead of a ForwardDiff result. This lets the outer ForwardDiff
# (around the steady-state solve, used as ground truth below) propagate
# without nested-Dual contamination.
F = AIBECSFunction((T_DIP, T_POP), (G_DIP, G_POP), nb, PmodelParameters; ∂Gs = ∂Gs)
@unpack DIP_geo = p
x = DIP_geo * ones(2nb) # initial guess

# Sanity-check the analytical ∂Gs against ForwardDiff before relying on
# them as the inner Jacobian. `check_∂Gs` returns a NamedTuple with
# `.maxabs` = max-abs difference; we expect it ≤ a few × eps on this
# small problem.
∂Gs_check = AIBECS.check_∂Gs((G_DIP, G_POP), ∂Gs, x, p, nb)
@test ∂Gs_check.maxabs < 1e-10

prob = SteadyStateProblem(F, x, p)
sol = solve(prob, CTKAlg())

# Observations from World Ocean Atlas to evaluate the AIBECS model
# mismatch with observations and generate the objective function.
ρSW = 1.035u"kg/L" # approximate mean sea water density
obs = let
    obs = WorldOceanAtlasTools.observations("phosphate")
    obs.value = ustrip.(upreferred.(obs.phosphate * ρSW))
    (obs,)
end
modify(DIP, POP) = (DIP,)
ωs = (1.0,) # weight for the DIP mismatch (POP weight = 0)
ωp = 1.0    # weight for the parameter prior penalty
f, ∇ₓf = f_and_∇ₓf(ωs, ωp, grd, modify, obs, PmodelParameters)

# Now apply the F-1 method.
# A tight abstol forces CTKAlg to re-converge after FiniteDiff's small
# parameter perturbations rather than returning the warm-started cached
# steady state.
solver_kwargs = (; abstol = 1e-12, maxItNewton = 200)
mem = F1Method.initialize_mem(F, ∇ₓf, x, λ, CTKAlg(); solver_kwargs...)
gradient(λ) = F1Method.gradient(f, F, ∇ₓf, mem, λ, CTKAlg(); solver_kwargs...)
hessian(λ)  = F1Method.hessian(f, F, ∇ₓf, mem, λ, CTKAlg(); solver_kwargs...)

# Cold-started reference objective for verification.
#
# F1Method's gradient/hessian use the analytical adjoint formula at the
# cached `mem.s`. To validate them, we need an independent ground truth
# that does *not* depend on the cached state. We re-solve from scratch
# (cold start) for each evaluation and let ForwardDiff propagate Duals
# through Newton + the analytical ∂Gs Jacobian. With analytical ∂Gs in
# AIBECS, this AD-through-the-solve is clean (no nested ForwardDiff).
function objective_cold(λ_eval)
    s_eval = solve(SteadyStateProblem(F, x, λ_eval), CTKAlg(); solver_kwargs...).u
    return f(s_eval, λ_eval)
end

# Mirror via the NonlinearSolve route. Independent IFT implementation
# (SciML's `__solve(::DualNonlinearProblem, …)` dispatch) — if both AD paths
# corroborate the analytical F1 formula, the formula is doubly cross-checked.
#
# `AIBECS.nonlinearproblem` attaches a *sparse* `jac_prototype`. SciML's
# nested-Dual (Hessian) path through `LinearSolveForwardDiffExt` currently
# fails for sparse Float64 prototype + Dual params (`MethodError:
# map!(partial_vals, ::Nothing, ::Vector{Float64})`). Dense / no prototype
# works. For this 2x2x2 test problem dense is fine — densifying lets us run
# both gradient *and* Hessian through the NL route end-to-end.
using SciMLBase: NonlinearProblem, NonlinearFunction
function objective_cold_nl(λ_eval)
    ssprob = SteadyStateProblem(F, x, λ_eval)
    nlprob = NonlinearProblem(ssprob)
    f_dense = NonlinearFunction(
        nlprob.f.f;
        jac = nlprob.f.jac,
        jac_prototype = Matrix{Float64}(undef, length(x), length(x)),
    )
    prob = NonlinearProblem(f_dense, nlprob.u0, nlprob.p)
    s_eval = solve(prob, NewtonRaphson(); solver_kwargs...).u
    return f(s_eval, λ_eval)
end

# F1 vs ForwardDiff + CTKAlg (existing checks)
@test ForwardDiff.gradient(objective_cold, 2λ) ≈ gradient(2λ) rtol = 1e-6
@test ForwardDiff.hessian(objective_cold, 2λ)  ≈ hessian(2λ)  rtol = 1e-6
@test ForwardDiff.gradient(objective_cold, λ)  ≈ gradient(λ)  rtol = 1e-6
@test ForwardDiff.hessian(objective_cold, λ)   ≈ hessian(λ)   rtol = 1e-6

# F1 vs ForwardDiff + NonlinearSolve
@test ForwardDiff.gradient(objective_cold_nl, 2λ) ≈ gradient(2λ) rtol = 1e-6
@test ForwardDiff.hessian(objective_cold_nl, 2λ)  ≈ hessian(2λ)  rtol = 1e-6
@test ForwardDiff.gradient(objective_cold_nl, λ)  ≈ gradient(λ)  rtol = 1e-6
@test ForwardDiff.hessian(objective_cold_nl, λ)   ≈ hessian(λ)   rtol = 1e-6

# FD reference (sanity check at looser rtol — second-order FD is noisy).
@test FiniteDiff.finite_difference_gradient(objective_cold, λ) ≈ gradient(λ) rtol = 1e-4
@test FiniteDiff.finite_difference_hessian(objective_cold, λ)  ≈ hessian(λ)  rtol = 1e-2

end # submodule
