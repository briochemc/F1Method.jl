 
mutable struct Mem 
    s     # 𝒔(𝒑)
    A     # factors of 𝐀 = ∇ₓ𝑭(𝒔,𝒑)
    ∇s    # ∇𝒔(𝒑)
    ∇ₓf   # ∇ₓ𝑓(𝒔,𝒑)
    p     # 𝒑
end

update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...) = if p ≠ mem.p 
    update_solution!(F, ∇ₓF, mem, p, alg; options...)
    s, m = mem.s.u, length(p)
    ∇ₚF = reduce(hcat, [𝔇(F(s, p + ε * e(j,m))) for j in 1:m]) # Eq.(18)
    mem.A = factorize(∇ₓF(s,p))   # update factors of ∇ₓ𝑭(𝒔,𝒑)
    mem.∇s .= mem.A \ -∇ₚF        # update ∇𝒔 via Eq.(13)
    mem.∇ₓf .= ∇ₓf(s,p)           # update ∇ₓ𝑓(𝒔,𝒑)
    mem.p = p                     # update 𝒑
end

function update_solution!(F, ∇ₓF, mem, p, alg; options...)
    if ~(mem.s isa SteadyStateSolution) || p ≠ mem.s.prob.p
        mem.s isa SteadyStateSolution ? x = mem.s.u : x = mem.s
        prob = SteadyStateProblem(F, ∇ₓF, x, p)    # define problem
        mem.s = solve(prob, alg; options...)       # update 𝒔
    end
end

function objective(f, F, ∇ₓF, mem, p, alg; options...)
    update_solution!(F, ∇ₓF, mem, p, alg; options...)
    return f(mem.s,p)
end

function gradient(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, ∇s, m = mem.s, mem.∇s, length(p)
    ∇ₚf = [𝔇(f(s,p + ε * e(j,m))) for j in 1:m]'    # Eq.(17)
    return mem.∇ₓf * ∇s + ∇ₚf                       # Eq.(12)
end

function hessian(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    update_mem!(f, F, ∇ₓf, ∇ₓF, mem, p, alg; options...)
    s, A, ∇s, m = mem.s, mem.A, mem.∇s, length(p)
    A⁻ᵀ∇ₓfᵀ = vec(A' \ mem.∇ₓf') # independent of (𝑗,𝑘)
    H = zeros(m,m)               # preallocate Hessian matrix
    for j in 1:m, k in j:m       # loop upper triangle (symmetry)
        pⱼₖ = p + ε₁ * e(j,m) + ε₂ * e(k,m)           # hyperdual 𝒑
        xⱼₖ = s + ε₁ * ∇s * e(j,m) + ε₂ * ∇s * e(k,m) # hyperdual 𝒙
        H[j,k] = ℌ(f(xⱼₖ,pⱼₖ)) - ℌ(F(xⱼₖ,pⱼₖ))' * A⁻ᵀ∇ₓfᵀ # Eq.(19)
        j ≠ k ? H[k,j] = H[j,k] : nothing # Hessian symmetry
    end
    return H
end

e(j, m) = [i == j for i in 1:m]      # 𝑗th basis vector of ℝᵐ
𝔇(x) = DualNumbers.dualpart.(x)      # dual part
ℌ(x) = HyperDualNumbers.ε₁ε₂part.(x) # hyperdual part
