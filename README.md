
<img src="https://user-images.githubusercontent.com/4486578/57202054-3d1c4400-6fe4-11e9-97d7-9a1ffbfcb2fc.png" alt="logo" title="F1method" align="right" height="200"/>

F-1 Method
==========

<p>
  <a href="https://doi.org/10.5281/zenodo.2667835">
    <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.2667835.svg" alt="DOI">
  </a>
  <a href="https://briochemc.github.io/F1Method.jl/stable/">
    <img src=https://img.shields.io/badge/docs-stable-blue.svg>
  </a>
</p>
<p>
  <a href="https://travis-ci.com/briochemc/F1Method.jl">
    <img alt="Build Status" src="https://travis-ci.com/briochemc/F1Method.jl.svg?branch=master">
  </a>
  <a href="https://ci.appveyor.com/project/briochemc/f1method-jl">
    <img alt="Build Status" src="https://ci.appveyor.com/api/projects/status/prm2xfd6q5pba1om?svg=true">
  </a>
</p>
<p>
  <a href='https://coveralls.io/github/briochemc/F1Method.jl'>
    <img src='https://coveralls.io/repos/github/briochemc/F1Method.jl/badge.svg' alt='Coverage Status' />
  </a>
  <a href="https://codecov.io/gh/briochemc/F1Method.jl">
    <img src="https://codecov.io/gh/briochemc/F1Method.jl/branch/master/graph/badge.svg" />
  </a>
</p>
<p>
  <a href="https://github.com/briochemc/F1Method.jl/blob/master/LICENSE">
    <img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-yellow.svg">
  </a>
</p>

This package implements the F-1 method descibed in Pasquier et al. (2019).
It allows for efficient quasi-auto-differentiation of an objective function defined implicitly by the solution of a steady-state problem represented by a discretized system of nonlinear partial differential equations that takes the form

<img src="https://latex.codecogs.com/svg.latex?&space;\boldsymbol{F}(\boldsymbol{x},&space;\boldsymbol{p})&space;=&space;0" title="Eq1"/>

where
<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{x}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{x}" title="\boldsymbol{x}" /></a>
is a column vector of the model state variables and
<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{p}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{p}" title="\boldsymbol{p}" /></a>
is a vector of parameters.
Specifically, the F-1 method allows for an efficient computation of both the gradient vector and the Hessian matrix of a generic objective function defined by
<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{f}(\boldsymbol{p})&space;\equiv&space;f(\boldsymbol{s},\boldsymbol{p})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{f}(\boldsymbol{p})&space;\equiv&space;f(\boldsymbol{s},\boldsymbol{p})" title="\hat{f}(\boldsymbol{p}) \equiv f(\boldsymbol{s},\boldsymbol{p})" /></a>
where
<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{s}(\boldsymbol{p})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{s}(\boldsymbol{p})" title="\boldsymbol{s}(\boldsymbol{p})" /></a>
is the steady-state solution of the system, i.e., such that
<a href="https://www.codecogs.com/eqnedit.php?latex=\boldsymbol{F}(\boldsymbol{s},&space;\boldsymbol{p})&space;=&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\boldsymbol{F}(\boldsymbol{s},&space;\boldsymbol{p})&space;=&space;0" title="\boldsymbol{F}(\boldsymbol{s}, \boldsymbol{p}) = 0" /></a>.

A requirement of the F-1 method is that the Jacobian matrix <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{A}&space;=&space;\nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x},\boldsymbol{p})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{A}&space;=&space;\nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x},\boldsymbol{p})" title="\mathbf{A} = \nabla_{\boldsymbol{x}}\boldsymbol{F}(\boldsymbol{x},\boldsymbol{p})" /></a> can be created, stored, and factorized.

To use the F-1 method, the user must:

- Make sure that there is a suitable algorithm `alg` to solve the steady-state equation
- overload the `solve` function and the `SteadyStateProblem` constructor from [DiffEqBase](https://github.com/JuliaDiffEq/DiffEqBase.jl). (An example is given in the CI tests — see, e.g., the [`test/simple_setup.jl`](test/simple_setup.jl) file.)
- Provide the derivatives of `f` and `F` with respect to the state, `x`.

### Simple usage

Make sure you have olverloaded `solve` from DiffEqBase.
Once an initial state, `x₀`, and some parameters, `p₀`, are chosen, simply evaluate the derivatives with

```julia
# Initialize the cache for storing reusable objects
mem = F1Method.initialize_mem(x₀, p₀)

# Wrap the objective, gradient, and Hessian functions
objective(p) = F1Method.f̂(f, F, ∇ₓF, mem, p, myAlg(); my_options...)
gradient(p) = F1Method.∇f̂(f, F, ∇ₓf, ∇ₓF, mem, p₀, myAlg(); my_options...)
hessian(p) = F1Method.∇²f̂(f, F, ∇ₓf, ∇ₓF, mem, p₀, myAlg(); my_options...)

# Compute the objective function, 𝑓̂(𝒑)
objective(p₀)

# Compute the gradient, ∇𝑓̂(𝒑)
gradient(p₀)

# Compute the Hessian matrix, ∇²𝑓̂(𝒑)
hessian(p₀)
```
