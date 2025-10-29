# `EpsteinLib.jl`

[![CI](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml)

This is a julia interface for [https://github.com/epsteinlib/epsteinlib](https://github.com/epsteinlib/epsteinlib) by Andreas A. Buchheit, Jonathan Busse, Ruben Gutendorf. Compiled binaries are available through [Epsteinlib_jll](https://github.com/JuliaBinaryWrappers/Epsteinlib_jll.jl).

The Epstein zeta function is defined by

```math
Z_{\nu, A}(x, y) = \sum_{\substack{z \in A \mathbb{Z}^d \\ z \ne x}} \frac{e^{-2\pi i y \cdot z}}{|x-z|^\nu},
```

if $\nu > d$, and the meromorphic continuation otherwise.

The `C` library may be called directly using the syntax
```julia
epsteinzeta(ν::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})::Complex{Float64}
```
and with optional keyword arguments as
```julia
epsteinzeta(ν; d, x, y, A)
```
where `d, x, y, A` are optional. `x` and `y` default to zero of size `d`, and `A` to the identity matrix of size `d`.  It performs some initial checks. 


The regularized Epstein zeta function is implemented in this library as
```julia
epsteinzetareg(ν::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})::Complex{Float64}
```
and with optional keyword arguments as
```julia
epsteinzetareg(ν; d, x, y, A)
```
