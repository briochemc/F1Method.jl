module FormulaOneMethod

using LinearAlgebra
using DualNumbers, HyperDualNumbers, DiffEqBase

mutable struct Mem
    s   # s(p)
    A   # factors of ∇ₓF(s, p)
    ∇s  # ∇s(p)
    ∇ₓf # ∇ₓf(s, p)
    p   # p
end

function update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    if p ≠ mem.p       # only update if p has changed
        update_solution!(F, ∇ₓF, mem, p, alg; options...)
        s, m = mem.s.u, length(p)
        ∇ₚF = hcat([𝔇(F(s, p + ε * e(j,m))) for j in 1:m]...) # Eq.(?)
        mem.A = factorize(∇ₓF(s,p))  # update factors of ∇ₓF(s, p)
        mem.∇s .= mem.A \ -∇ₚF    # update ∇s via Eq.(?)
        mem.∇ₓf .= ∇ₓf(s,p)          # update ∇ₓf(s, p)
        mem.p = p      # update p
    end
end

function update_solution!(F, ∇ₓF, mem, p, alg; options...)
    if ~(mem.s isa SteadyStateSolution) || p ≠ mem.s.prob.p
        mem.s isa SteadyStateSolution ? x = mem.s.u : x = mem.s
        prob = SteadyStateProblem(F, ∇ₓF, x, p) # define problem
        mem.s = solve(prob, alg; options...) # update s (inner solver)
    end
end

function f̂(f, F, ∇ₓF, mem, p, alg; options...) # objective
    update_solution_only!(F, ∇ₓF, mem, p, alg; options...)
    return f(mem.s,p)
end

function ∇f̂(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...) # gradient
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, ∇s, m = mem.s, mem.∇s, length(p)
    ∇ₚf = [𝔇(f(s,p + ε * e(j,m))) for j in 1:m]' # Eq. (?)
    return mem.∇ₓf * ∇s + ∇ₚf               # Eq. (?)
end

function ∇²f̂(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...) # Hessian
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, A, ∇s, m = mem.s, mem.A, mem.∇s, length(p)
    A⁻ᵀ∇ₓfᵀ = vec(A' \ mem.∇ₓf') # independent of (j,k)
    out = zeros(m,m)       # preallocate
    for j in 1:m, k in j:m # Loop for Eq.(?)
        pⱼₖ = p + ε₁ * e(j,m) + ε₂ * e(k,m)           # Hyperdual p
        xⱼₖ = s + ε₁ * ∇s * e(j,m) + ε₂ * ∇s * e(k,m) # Hyperdual x
        out[j,k] = ℌ(f(xⱼₖ,pⱼₖ)) - ℌ(F(xⱼₖ,pⱼₖ))' * A⁻ᵀ∇ₓfᵀ # Eq.(?)
        j ≠ k ? out[k,j] = out[j,k] : nothing   # symmetry
    end
    return out
end

function initialize_mem(x, p)
    n, m = length(x), length(p)
    return Mem(copy(x), nothing, zeros(n,m), zeros(1,n), nothing)
end

# Helper functions
e(j, m) = [i == j for i in 1:m]      # j-th basis vector
𝔇(x) = DualNumbers.dualpart.(x)      # dual part
ℌ(x) = HyperDualNumbers.ε₁ε₂part.(x) # hyperdual part

end # module
