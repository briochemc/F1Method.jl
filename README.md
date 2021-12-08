
<img src="https://user-images.githubusercontent.com/4486578/57202054-3d1c4400-6fe4-11e9-97d7-9a1ffbfcb2fc.png" alt="logo" title="F1method" align="right" height="200"/>

# F-1 algorithm

<p>
  <a href="https://briochemc.github.io/F1Method.jl/stable/">
    <img src="https://img.shields.io/github/workflow/status/briochemc/F1Method.jl/Documentation?style=for-the-badge&label=Documentation&logo=Read%20the%20Docs&logoColor=white">
  </a>
</p>

<p>
  <a href="https://doi.org/10.5281/zenodo.2667835">
    <img src="http://img.shields.io/badge/DOI-10.5281%20%2F%20zenodo.2667835-blue.svg?&style=flat-square">
  </a>
  <a href="https://github.com/briochemc/F1Method.jl/blob/master/LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-blue.svg?&style=flat-square">
  </a>
</p>

<p>
  <a href="https://github.com/briochemc/F1Method.jl/actions">
    <img src="https://img.shields.io/github/workflow/status/briochemc/F1Method.jl/Mac%20OS%20X?label=OSX&logo=Apple&logoColor=white&style=flat-square">
  </a>
  <a href="https://github.com/briochemc/F1Method.jl/actions">
    <img src="https://img.shields.io/github/workflow/status/briochemc/F1Method.jl/Linux?label=Linux&logo=Linux&logoColor=white&style=flat-square">
  </a>
  <a href="https://github.com/briochemc/F1Method.jl/actions">
    <img src="https://img.shields.io/github/workflow/status/briochemc/F1Method.jl/Windows?label=Windows&logo=Windows&logoColor=white&style=flat-square">
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
- overload the `solve` function and the `SteadyStateProblem` constructor from [SciMLBase](https://github.com/JuliaDiffEq/SciMLBase.jl). (An example is given in the CI tests — see, e.g., the [`test/simple_setup.jl`](test/simple_setup.jl) file.)
- Provide the derivatives of `f` and `F` with respect to the state, `x`.

## A concrete example

Make sure you have overloaded `solve` from SciMLBase
(an example of how to do this is given in the [documentation](https://briochemc.github.io/F1Method.jl/stable/)).
Once initial values for the state, `x`, and parameters, `p`, are chosen, simply initialize the required memory cache, `mem` via

```julia
# Initialize the cache for storing reusable objects
mem = initialize_mem(F, ∇ₓf, x, p, alg; options...)
```

wrap the functions into functions of `p` only via

```julia
# Wrap the objective, gradient, and Hessian functions
objective(p) = F1Method.objective(f, F, mem, p, alg; options...)
gradient(p) = F1Method.gradient(f, F, ∇ₓf, mem, p, alg; options...)
hessian(p) = F1Method.hessian(f, F, ∇ₓf, mem, p, alg; options...)
```

and compute the objective, gradient, or Hessian via either of

```julia
objective(p)

gradient(p)

hessian(p)
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