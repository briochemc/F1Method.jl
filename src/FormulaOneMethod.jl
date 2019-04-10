module FormulaOneMethod

using LinearAlgebra, DualNumbers, HyperDualNumbers, DiffEqBase
 
mutable struct Mem # Storage for efficient reuse
    s     # 𝑠(𝑝)
    A     # factors of A = ∇ₓ𝐹(𝑠,𝑝)
    ∇s    # ∇𝑠(𝑝)
    ∇ₓf   # ∇ₓ𝑓(𝑠,𝑝)
    p     # 𝑝
end

function update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    if p ≠ mem.p                    # only update mem if 𝑝 has changed
        update_solution!(F, ∇ₓF, mem, p, alg; options...)
        s, m = mem.s.u, length(p)
        ∇ₚF = hcat([𝔇(F(s, p + ε * e(j,m))) for j in 1:m]...) # Eq.(?)
        mem.A = factorize(∇ₓF(s,p))   # update factors of ∇ₓ𝐹(𝑠,𝑝)
        mem.∇s .= mem.A \ -∇ₚF        # update ∇𝑠               Eq.(?)
        mem.∇ₓf .= ∇ₓf(s,p)           # update ∇ₓ𝑓(𝑠,𝑝)
        mem.p = p                     # update 𝑝
    end
end

function update_solution!(F, ∇ₓF, mem, p, alg; options...)
    if ~(mem.s isa SteadyStateSolution) || p ≠ mem.s.prob.p
        mem.s isa SteadyStateSolution ? x = mem.s.u : x = mem.s
        prob = SteadyStateProblem(F, ∇ₓF, x, p)       # define problem
        mem.s = solve(prob, alg; options...)          # update 𝑠
    end
end

function f̂(f, F, ∇ₓF, mem, p, alg; options...)             # objective
    update_solution_only!(F, ∇ₓF, mem, p, alg; options...)
    return f(mem.s,p)
end

function ∇f̂(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)        # gradient
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, ∇s, m = mem.s, mem.∇s, length(p)
    ∇ₚf = [𝔇(f(s,p + ε * e(j,m))) for j in 1:m]'    # Eq.(?)
    return mem.∇ₓf * ∇s + ∇ₚf                       # Eq.(?)
end

function ∇²f̂(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)        # Hessian
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, A, ∇s, m = mem.s, mem.A, mem.∇s, length(p)
    A⁻ᵀ∇ₓfᵀ = vec(A' \ mem.∇ₓf') # independent of (𝑗,𝑘)
    H = zeros(m,m)               # preallocate Hessian matrix
    for j in 1:m, k in j:m       # loop upper triangle (symmetry)
        pⱼₖ = p + ε₁ * e(j,m) + ε₂ * e(k,m)              # hyperdual 𝑝
        xⱼₖ = s + ε₁ * ∇s * e(j,m) + ε₂ * ∇s * e(k,m)    # hyperdual 𝑥
        H[j,k] = ℌ(f(xⱼₖ,pⱼₖ)) - ℌ(F(xⱼₖ,pⱼₖ))' * A⁻ᵀ∇ₓfᵀ     # Eq.(?)
        j ≠ k ? H[k,j] = H[j,k] : nothing
    end
    return H
end

e(j, m) = [i == j for i in 1:m]      # 𝑗ᵗʰ basis vector of ℝᵐ
𝔇(x) = DualNumbers.dualpart.(x)      # dual part
ℌ(x) = HyperDualNumbers.ε₁ε₂part.(x) # hyperdual part

function initialize_mem(x, p)             # function to initialize mem
    n, m = length(x), length(p)
    return Mem(copy(x), nothing, zeros(n,m), zeros(1,n), nothing)
end

end
