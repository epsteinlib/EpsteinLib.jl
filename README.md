# `EpsteinLib.jl`
<img align="right" src="https://avatars.githubusercontent.com/u/177750891?v=4" width="110">

[![CI](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml)


Authors: David Gómez-Castro and Jonathan K. Busse

Julia interface for the C library [epsteinlib](https://github.com/epsteinlib/epsteinlib) by Andreas A. Buchheit, Jonathan K. Busse, and Ruben Gutendorf.

Precompiled binaries are available through [Epsteinlib_jll](https://github.com/JuliaBinaryWrappers/Epsteinlib_jll.jl).

For a $d$-dimensional lattice $\Lambda=A\mathbb Z^d$, with $A\in \mathbb R^{d\times d}$ regular, $\boldsymbol x,\boldsymbol y \in \mathbb R^d$, and $\nu \in \mathbb C$, the Epstein zeta function is defined by the Dirichlet series

$$
Z_{\Lambda,\nu}(\boldsymbol x,\boldsymbol y)
= \sum_{z \in \Lambda}{}^{'} \frac{e^{-2\pi i \boldsymbol y \cdot \boldsymbol z}}{\left| \boldsymbol x- \boldsymbol z\right|^\nu},\quad \mathrm{Re}(\nu)>d,
$$

which can be meromorphically continued to $\nu \in \mathbb C$. Here, the primed sum excludes the case $\boldsymbol z = \boldsymbol x.$

The Epstein zeta function is implemented as
```julia
epsteinzeta(ν::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})::Complex{Float64}
```
and with optional keyword arguments as
```julia
epsteinzeta(ν; d, A, x, y)
```
where at least one of the arguments `d`, `x`, `y`, or `A` must be provided. By default, `x` and `y` are zero vectors of length `d`, and `A` is the `d × d` identity matrix.

In addition, this library includes the regularized Epstein zeta function, which is analytic around $\boldsymbol y= \boldsymbol 0$, and is defined via

$$
Z_{\Lambda,\nu}^{\mathrm{reg}}(\boldsymbol x,\boldsymbol y) =
e^{2\pi i \boldsymbol x\cdot\boldsymbol y}
Z_{\Lambda,\nu}(\boldsymbol x,\boldsymbol y )
-\frac{\hat{s}_{\nu}(\boldsymbol y)}{V_{\Lambda}},
$$

where $V_{\Lambda}=|\det A|$ is the volume of the elementary lattice cell, and the Fourier transform of the singularity $s_{\nu}=|\boldsymbol{\cdot}|^{-\nu}$ is defined as in [epsteinlib](https://github.com/epsteinlib/epsteinlib).

The regularized Epstein zeta function is implemented as
```julia
epsteinzetareg(ν::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})::Complex{Float64}
```
and with optional keyword arguments as
```julia
epsteinzetareg(ν; d, A, x, y)
```
Defaults for `x`, `y`, and `A` are identical to those used in `epsteinzeta`.

## Anisotropic Epstein zeta function

Let $\nu\in\mathbb C$ and signify by the multi-index $\boldsymbol\alpha\in\mathbb N_0^d$ the anisotropy strength of

$$
V_{\nu,\boldsymbol \alpha}(\boldsymbol z)
= \frac{\boldsymbol z^{\boldsymbol \alpha}}{\vert \boldsymbol z \vert^\nu}
,\qquad
\boldsymbol z\in\mathbb R^d\setminus\{\boldsymbol 0\},
$$

with $\boldsymbol z^{\boldsymbol\alpha}=z_1^{\alpha_1}z_2^{\alpha_2}\ldots z_d^{\alpha_d}$. For a $d$-dimensional lattice $\Lambda$ and $\boldsymbol x,\boldsymbol y \in \mathbb R^d$, the anisotropic Epstein zeta function is then define as

$$
Z_{\Lambda,\nu,\boldsymbol \alpha}(\boldsymbol x,\boldsymbol y)
= \sum_{z \in \Lambda}{}^{'} e^{-2\pi i \boldsymbol y \cdot \boldsymbol z}V_{\nu,\boldsymbol \alpha}(\boldsymbol z-\boldsymbol x),\quad \mathrm{Re}(\nu)>d
+|\boldsymbol \alpha|,
$$

meromorphically continued to $\nu \in \mathbb C$; where we define $|\boldsymbol{\alpha}|=\alpha_1+\ldots+\alpha_d$.  Here, we recover the Epstein zeta function for $\boldsymbol \alpha=\boldsymbol 0$. 

The anisotropic Epstein zeta function is implemented as

```julia
epsteinzetaaniso(ν::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64}, α::Vector{UInt32})::Complex{Float64}
```

and with optional keyword arguments as

```julia
epsteinzetaaniso(ν, α; d, A, x, y)
```

where `α` determines the dimension if none of `d`, `x`, `y`, or `A` is provided. Defaults for `x`, `y`, and `A` are identical to those used in `epsteinzeta`.

In addition, the library includes the regularized anisotropic Epstein zeta function defined via

$$
Z_{\Lambda, \nu,\boldsymbol\alpha}^{(\mathrm{reg})}(\boldsymbol x,\boldsymbol y) = e^{2\pi i \boldsymbol{x}\cdot\boldsymbol{y}}Z_{\Lambda,\nu,\boldsymbol\alpha}(\boldsymbol x,\boldsymbol y) -\frac{\hat s^{(\boldsymbol\alpha)}_{\nu}(\boldsymbol y)}{(-2\pi i)^{|\boldsymbol\alpha|}V_{\Lambda}}
,\qquad \boldsymbol y\neq \boldsymbol 0,
$$

and continuously extended to $\boldsymbol y=\boldsymbol 0$, where $`\hat s^{(\boldsymbol\alpha)}_\nu`$ denotes the $\boldsymbol\alpha$-derivative of $\hat{s}_\nu$.

The regularized anisotropic Epstein zeta function is implemented as

```julia
epsteinzetaanisoreg(ν::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64}, α::Vector{UInt32})::Complex{Float64}
```

and with optional keyword arguments as

```julia
epsteinzetaanisoreg(ν, α; d, A, x, y)
```

Defaults for `x`, `y`, and `A` are identical to those used in `epsteinzeta`.

## Installation and usage

The library can be installed via

```julia
using Pkg; Pkg.add("EpsteinLib") # Stable release (recommended)

# Development / latest GitHub version (optional)
# Pkg.add(url="https://github.com/epsteinlib/EpsteinLib.jl")
```

The following example then computes the Madelung constant to machine precision 

```julia
using EpsteinLib, Printf

ν = 1.0
A = [1.0 0.0 0.0;
     0.0 1.0 0.0;
     0.0 0.0 1.0]
x = [0.0, 0.0, 0.0]
y = [0.5, 0.5, 0.5]

# Calculate Madelung constant
madelung = epsteinzeta(ν, A, x, y)

# Reference value and relative error
madelung_ref = -1.7475645946331821906362120355443974
relerr = abs(madelung_ref - real(madelung)) / abs(madelung_ref)

println("Madelung sum in 3 dimensions:\t", real(madelung))
println("Reference value:\t\t", madelung_ref)
@printf("Relative error:\t\t\t +%.2e\n", relerr)
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please use [pre-commit](https://pre-commit.com) to ensure your commits are well-formatted by running 
```sh
pip install pre-commit
pre-commit install
```
when you clone the repo.
