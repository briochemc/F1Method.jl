
<img src="https://user-images.githubusercontent.com/4486578/57202054-3d1c4400-6fe4-11e9-97d7-9a1ffbfcb2fc.png" alt="logo" title="F1method" align="right" height="200"/>

# F-1 algorithm

<p>
  <a href="https://doi.org/10.5281/zenodo.2667835">
    <img src="http://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.2667835-blue.svg?&style=flat-square">
  </a>
  <a href="https://github.com/briochemc/F1Method.jl/blob/master/LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-blue.svg?&style=flat-square">
  </a>
</p>

<p>
  <a href="https://github.com/briochemc/F1Method.jl/actions/workflows/Test.yml">
    <img src="https://img.shields.io/github/actions/workflow/status/briochemc/F1Method.jl/Test.yml?label=Test&logo=github&logoColor=white&style=flat-square">
  </a>
  <a href="https://codecov.io/gh/briochemc/F1Method.jl">
    <img src="https://img.shields.io/codecov/c/github/briochemc/F1Method.jl/master?label=Codecov&logo=codecov&logoColor=white&style=flat-square">
  </a>
</p>

This package implements the F-1 algorithm described in *[Pasquier and Primeau](https://www.bpasquier.com/publication/pasquier_primeau_sisc_2019/)* (in preparation).
It allows for efficient quasi-auto-differentiation of an objective function defined implicitly by the solution of a steady-state problem.

Consider a discretized system of nonlinear partial differential equations that takes the form

```
F(x,p) = 0
```

where `x` is a column vector of the model state variables and `p` is a vector of parameters.
The F-1 algorithm then allows for an efficient computation of both the gradient vector and the Hessian matrix of a generic objective function defined by

```
objective(p) = f(s(p),p)
```

where `s(p)` is the steady-state solution of the system, i.e., such that `F(s(p),p) = 0` and where `f(x,p)` is for example a measure of the mismatch between observed state, parameters, and observations.
Optimizing the model is then simply done by minimizing `objective(p)`.
(See *[Pasquier and Primeau](https://www.bpasquier.com/publication/pasquier_primeau_sisc_2019/)* (in preparation), for more details.)

## Advantages of the F-1 algorithm

The F-1 algorithm is **easy** to use, gives **accurate** results, and is computationally **fast**:

- **Easy** — The F-1 algorithm basically just needs the user to provide a solver (for finding the steady-state), the mismatch function, `f`, an ODEFunction, `F` with its Jacobian, and the gradient of the objective w.r.t. `∇ₓf`.
    (Note these derivatives can be computed numerically, via the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package for example.)
- **Accurate** — Thanks to ForwardDiff's nested dual numbers implementation, the accuracy of the gradient and Hessian, as computed by the F-1 algorithm, are close to machine precision.
- **Fast** — The F-1 algorithm is as fast as if you derived analytical formulas for every first and second derivatives *and* used those in the most efficient way.
    This is because the bottleneck of such computations is the number of matrix factorizations, and the F-1 algorithm only requires a single one. In comparison, standard autodifferentiation methods that take the steady-state solver as a black box would require order `m` or `m^2` factorizations, where `m` is the number of parameters.

## What's needed?

A requirement of the F-1 algorithm is that the Jacobian matrix `A = ∇ₓF` can be created, stored, and factorized.

To use the F-1 algorithm, the user must:

- Make sure that there is a suitable algorithm `alg` to solve the steady-state equation
- overload the `solve` function and the `SteadyStateProblem` constructor from [SciMLBase](https://github.com/SciML/SciMLBase.jl). (An example is given below, and another in the CI tests — see, e.g., the [`test/simple_setup.jl`](test/simple_setup.jl) file.)
- Provide the derivatives of `f` and `F` with respect to the state, `x`.

## A concrete example

We define a state function, `F(x,p)`, to which we apply a Newton solver to find the steady-state solution, `x`, such that `F(x,p) = 0`.
This defines the steady-state solution as an implicit function of the parameters, `p`, which we denote by `s(p)`.
The Newton solver requires the Jacobian, `∇ₓF`, to update the state iterates, so we start by creating the functions `F(x,p)` and `∇ₓF(x,p)`.
As an example, we use a simple model with two state variables and two parameters.
(For simplicity we use [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) to evaluate the Jacobian.)

```julia
using F1Method
using LinearAlgebra, SciMLBase, ForwardDiff

# State function F
# (The `t = 0` default lets us call `statefun(x, p)` directly while still
# satisfying SciMLBase's ODE-shape signature check at `ODEFunction` construction.
# `ODEFunction{false}` declares it out-of-place and skips the inplace probe.)
statefun(x, p, t = 0) = [
    -2 * (p[1] - x[1]) - 4 * p[2] * (x[2] - x[1]^2) * x[1]
    p[2] * (x[2] - x[1]^2)
]
F = ODEFunction{false}(statefun;
    jac = (x, p, t = 0) -> ForwardDiff.jacobian(x -> statefun(x, p), x))
```

We also define a cost function `f(x,p)` (that we wish to minimize under the constraint that `F(x,p) = 0`).
The F-1 method requires the derivative of `f` w.r.t. the state, `x`, so we use ForwardDiff again:

```julia
# Define mismatch function f(x,p) and its derivative ∇ₓf(x,p)
function state_mismatch(x)
    δ(x) = x .- 1
    return 0.5δ(x)'δ(x)
end
function parameter_mismatch(p)
    δ(p) = log.(p)
    return 0.5δ(p)'δ(p)
end
f(x,p) = state_mismatch(x) + parameter_mismatch(p)
∇ₓf(x,p) = ForwardDiff.jacobian(x -> [f(x,p)], x)
```

Once these are set up, we tell the F-1 method how to solve for the steady state by using the [SciMLBase](https://github.com/SciML/SciMLBase.jl) API: we write a small Newton solver, then overload `solve` and the `SteadyStateProblem` constructor.

```julia
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
    SciMLBase.build_solution(prob, alg, x_steady, resid; retcode=:Success)
end
```

Pick initial values for the state `x` and parameters `p`, initialize the cache, and wrap the objective, gradient, and Hessian:

```julia
x₀, p₀ = [1.0, 2.0], [3.0, 4.0]

# Initialize the cache for storing reusable objects
mem = F1Method.initialize_mem(F, ∇ₓf, x₀, p₀, MyAlg())

# Define the functions via the F1 method
objective(p) = F1Method.objective(f, F, mem, p, MyAlg())
gradient(p)  = F1Method.gradient(f, F, ∇ₓf, mem, p, MyAlg())
hessian(p)   = F1Method.hessian(f, F, ∇ₓf, mem, p, MyAlg())
```

You can now call `objective(p₀)`, `gradient(p₀)`, and `hessian(p₀)` directly:

```julia
julia> objective(p₀)
35.56438050824269

julia> gradient(p₀)
1×2 Array{Float64,2}:
 50.3662  0.346574

julia> hessian(p₀)
2×2 Array{Float64,2}:
 52.989   0.0
  0.0    -0.0241434
```

That's it.
You were told it was simple, weren't you?
Now you can test how fast and accurate it is!

## Citing the software

If you use this package, or implement your own package based on the F-1 algorithm please cite us.
If you use the F-1 algorithm, please cite *[Pasquier and Primeau](https://www.bpasquier.com/publication/pasquier_primeau_sisc_2019/)* (in prep.).
If you also use this package directly, please cite it! (Use [the Zenodo link](https://doi.org/10.5281/zenodo.2667835) or the [CITATION.bib file](./CITATION.bib), which contains a bibtex entry.)

# Future

This package is developed mainly for use with [AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl) and is likely not in its final form.
The API was just changed in v0.5 (to match the API changes in AIBECS.jl v0.11).
That being said, ultimately, it would make sense for the shortcuts used here to be integrated into a package like ChainRules.jl.
For the time being, AIBECS users can use F1Method.jl to speed up their optimizations.
