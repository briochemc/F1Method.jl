module F1MethodLinearSolveExt

using F1Method
using F1Method: F1Cache
using LinearSolve
using SciMLBase
using DifferentiationInterface
const DI = DifferentiationInterface

function F1Method._init_F1Cache(linsolve::SciMLBase.AbstractLinearAlgorithm,
                                F, ∇ₓf, x, p, alg, ad; options...)
    x = copy(x)
    p = copy(p)
    psol = copy(p)
    s = solve(SteadyStateProblem(F, x, p), alg; options...).u
    J = F.jac(s, p)
    ∇ₚF = DI.jacobian(p -> F(s, p), ad, p)
    # `u0 = similar(∇ₚF)` forces LinearCache.u to match the matrix-RHS
    # shape; without it LinearSolve.init allocates `u` as a Vector even
    # when `b` is a Matrix, and the subsequent solve! errors on the
    # shape mismatch.
    lc = init(LinearProblem(J, -∇ₚF; u0 = similar(∇ₚF)), linsolve)
    sol = solve!(lc)
    ∇s = copy(sol.u)                         # F1Cache owns its own ∇s buffer
    return F1Cache(lc, s, ∇s, ∇ₓf(s, p), p, psol, ad)
end

# cache.linear_cache is a LinearSolve.LinearCache; .A is the cache's
# stored Jacobian. Assigning to it triggers LinearSolve's setproperty!
# → isfresh = true, so the next solve!(linear_cache) refactors.
# UMFPACK/KLU reuse the symbolic factorization automatically when
# colptr/rowval match — that's the main win beyond plain `\`.
F1Method.set_A!(cache::F1Cache{<:LinearSolve.LinearCache}, A) =
    cache.linear_cache.A = A

# solve!(lc).u aliases lc.u. The only call site is
# `cache.∇s .= solveA!(...)`, which copies — don't "optimize" the
# broadcast away or you'll alias cache.∇s to the LinearCache buffer.
function F1Method.solveA!(cache::F1Cache{<:LinearSolve.LinearCache}, b)
    cache.linear_cache.b = b
    return solve!(cache.linear_cache).u
end

# `cacheval` is the algorithm's stored Factorization (or, for
# GenericLUFactorization, a `(LU, ipiv)` tuple whose first element is
# the factorization). Group A algorithms make `factor' \ B` legal.
# Upstream-recommended idiom for concrete matrices until
# SciML/LinearSolve.jl#92 ships a first-class adjoint toggle on
# LinearCache.
_factor(cv) = cv
_factor(cv::Tuple) = first(cv)
F1Method.solveAᵀ!(cache::F1Cache{<:LinearSolve.LinearCache}, B) =
    _factor(cache.linear_cache.cacheval)' \ B

end # module
