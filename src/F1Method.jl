module F1Method

#======================================================================
This package (the code below) implements the F-1 method as described
in the work of Pasquier et al. (2019). The numbers in parentheses
refer to the Equation numbers in the above manuscript. A bibtex
citation file is available in the GitHub repository.
======================================================================#

using LinearAlgebra, DualNumbers, HyperDualNumbers, DiffEqBase

mutable struct Mem # memory cache storing reusable objects
    s     # 𝒔(𝒑)
    A     # factors of 𝐀 = ∇ₓ𝑭(𝒔,𝒑)
    ∇s    # ∇𝒔(𝒑)
    ∇ₓf   # ∇ₓ𝑓(𝒔,𝒑)
    p     # 𝒑
end

function update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    if p ≠ mem.p                      # only update mem if 𝒑 has changed
        update_solution!(F, ∇ₓF, mem, p, alg; options...)
        s, m = mem.s.u, length(p)
        ∇ₚF = reduce(hcat, [𝔇(F(s, p + ε * e(j,m))) for j in 1:m]) #(18)
        mem.A = factorize(∇ₓF(s,p))   # update factors of ∇ₓ𝑭(𝒔,𝒑)
        mem.∇s .= mem.A \ -∇ₚF        # update ∇𝒔 via (13)
        mem.∇ₓf .= ∇ₓf(s,p)           # update ∇ₓ𝑓(𝒔,𝒑)
        mem.p = p                     # update 𝒑
    end
end

function update_solution!(F, ∇ₓF, mem, p, alg; options...)
    if ~(mem.s isa SteadyStateSolution) || p ≠ mem.s.prob.p
        mem.s isa SteadyStateSolution ? x = mem.s.u : x = mem.s
        prob = SteadyStateProblem(F, ∇ₓF, x, p)       # define problem
        mem.s = solve(prob, alg; options...)          # update 𝒔
    end
end

function f̂(f, F, ∇ₓF, mem, p, alg; options...)             # objective
    update_solution!(F, ∇ₓF, mem, p, alg; options...)
    return f(mem.s,p)
end

function ∇f̂(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)        # gradient
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, ∇s, m = mem.s, mem.∇s, length(p)
    ∇ₚf = [𝔇(f(s,p + ε * e(j,m))) for j in 1:m]'    # (17)
    return mem.∇ₓf * ∇s + ∇ₚf                       # (12)
end

function ∇²f̂(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)        # Hessian
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, A, ∇s, m = mem.s, mem.A, mem.∇s, length(p)
    A⁻ᵀ∇ₓfᵀ = vec(A' \ mem.∇ₓf') # independent of (𝑗,𝑘)
    H = zeros(m,m)               # preallocate Hessian matrix
    for j in 1:m, k in j:m       # loop upper triangle (symmetry)
        pⱼₖ = p + ε₁ * e(j,m) + ε₂ * e(k,m)              # hyperdual 𝒑
        xⱼₖ = s + ε₁ * ∇s * e(j,m) + ε₂ * ∇s * e(k,m)    # hyperdual 𝒙
        H[j,k] = ℌ(f(xⱼₖ,pⱼₖ)) - ℌ(F(xⱼₖ,pⱼₖ))' * A⁻ᵀ∇ₓfᵀ    # (19)
        j ≠ k ? H[k,j] = H[j,k] : nothing
    end
    return H
end

e(j, m) = [i == j for i in 1:m]      # 𝑗th basis vector of ℝᵐ
𝔇(x) = DualNumbers.dualpart.(x)      # dual part
ℌ(x) = HyperDualNumbers.ε₁ε₂part.(x) # hyperdual part

function initialize_mem(x, p) # function to initialize the cache (mem)
    n, m = length(x), length(p)
    return Mem(copy(x), nothing, zeros(n,m), zeros(1,n), nothing)
end

end
