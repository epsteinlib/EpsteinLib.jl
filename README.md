# `EpsteinZetaFunction.jl`

[![CI](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml)

This is a julia interface for [https://github.com/epsteinlib/epsteinlib](https://github.com/epsteinlib/epsteinlib) by Andreas A. Buchheit, Jonathan Busse, Ruben Gutendorf.


The Epstein zeta function is implemented in this library as
```julia
epsteinzeta(ν::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})::Complex{Float64}
```

and with optional arguments as
```julia
epsteinzeta(
    ν::T0;
    d::Union{Integer,Nothing}=nothing,
    x::Union{Vector{T1},Nothing}=nothing,
    y::Union{Vector{T2},Nothing}=nothing,
    A::Union{Matrix{T3},Nothing}=nothing,
)::Complex{Float64} where {T0<:Real,T1<:Real,T2<:Real,T3<:Real}
```

The regularized Epstein zeta function is implemented in this library as
```julia
epsteinzetareg(ν::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})::Complex{Float64}
```

and with optional arguments as
```julia
epsteinzetareg(
    ν::T0;
    d::Union{Integer,Nothing}=nothing,
    x::Union{Vector{T1},Nothing}=nothing,
    y::Union{Vector{T2},Nothing}=nothing,
    A::Union{Matrix{T3},Nothing}=nothing,
)::Complex{Float64} where {T0<:Real,T1<:Real,T2<:Real,T3<:Real}

```