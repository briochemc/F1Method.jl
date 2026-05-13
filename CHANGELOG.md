# Changelog

## 0.7.0

Breaking renames (deprecation wrappers retained for 0.7.x; remove in 0.8):

- `Mem` → `F1Cache`
- `initialize_mem(...)` → `F1Cache(...)` constructor
- `update_mem!` → `update_cache!` (internal)
- `Mem.A` field → `F1Cache.linear_cache`

Migration:

```julia
# 0.6.x
mem = F1Method.initialize_mem(F, ∇ₓf, x₀, p₀, alg)
F1Method.objective(f, F, mem, p, alg)
F1Method.gradient(f, F, ∇ₓf, mem, p, alg)
F1Method.hessian(f, F, ∇ₓf, mem, p, alg)

# 0.7.0
cache = F1Method.F1Cache(F, ∇ₓf, x₀, p₀, alg)
F1Method.objective(f, F, cache, p, alg)
F1Method.gradient(f, F, ∇ₓf, cache, p, alg)
F1Method.hessian(f, F, ∇ₓf, cache, p, alg)
```

`initialize_mem(...)` still works in 0.7.x with a deprecation warning;
it will be removed in 0.8.

New:

- Optional LinearSolve.jl backend, exposed via the new `linsolve`
  kwarg on the `F1Cache` constructor. With it,
  `F1Cache.linear_cache` holds a `LinearSolve.LinearCache`, enabling
  symbolic factorization reuse for sparse direct solvers and buffer
  reuse for dense. Default behaviour (`linsolve === nothing`) is
  unchanged from 0.6.x.
- Tested with `LUFactorization`, `GenericLUFactorization`,
  `UMFPACKFactorization`, `KLUFactorization`.
