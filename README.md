# `EpsteinZetaFunction.jl`

[![CI](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml)

This is a julia interface for [https://github.com/epsteinlib/epsteinlib](https://github.com/epsteinlib/epsteinlib) by Andreas A. Buchheit, Jonathan Busse, Ruben Gutendorf.

The Epstein zeta function is defined by
$$
    Z_{\nu, A}(x, y) = \sum_{\substack{z \in A \mathbb{Z}^d \\ z \ne x}} \frac{-e^{2\pi i y \cdot z}}{|x-z|^\nu}.
$$
if the real part of nu is greater than the system dimension, and the meromorphic continuation otherwise.

It provides the convenient syntax
```julia
    epsteinzeta(ν; d, x, y, A)
```
where `d, x, y, A` are optional. `x` and `y` default to zero of size `d`, and `A` to the identity matrix of size `d`.  It performs some initial checks. 
The `C` library may be called directly using the syntax
```julia
    epsteinzeta(
        ν::Float64,
        A::Matrix{Float64},
        x::Vector{Float64},
        y::Vector{Float64},
    )::Complex{Float64}
```