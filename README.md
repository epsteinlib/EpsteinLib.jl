# `EpsteinLib.jl`
<img align="right" src="https://avatars.githubusercontent.com/u/177750891?v=4" width="110">

[![CI](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/dgomezcastro/EpsteinZetaFunction.jl/actions/workflows/ci.yml)


Authors: David Gómez-Castro and Jonathan K. Busse

Julia interface for the C library [epsteinlib](https://github.com/epsteinlib/epsteinlib) by Andreas A. Buchheit, Jonathan K. Busse, and Ruben Gutendorf.

Precompiled binaries are available through [Epsteinlib_jll](https://github.com/JuliaBinaryWrappers/Epsteinlib_jll.jl).

For a $d$-dimensional lattice $\Lambda=A\mathbb Z^d$, with $A\in \mathbb R^{d\times d}$ regular, $\boldsymbol x,\boldsymbol y \in \mathbb R^d$, and $\nu \in \mathbb C$, the Epstein zeta function is defined by the Dirichlet series

$$
Z_{\Lambda,\nu}\begin{vmatrix} \boldsymbol x \newline\boldsymbol y \end{vmatrix}
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

In addition, this wrapper includes the regularized Epstein zeta function, which is analytic around $\boldsymbol y=0$, and is defined via

$$
Z_{\Lambda,\nu}^{\mathrm{reg}}\begin{vmatrix} \boldsymbol x \newline\boldsymbol y \end{vmatrix} =
e^{2\pi i \boldsymbol x\cdot\boldsymbol y}
Z_{\Lambda,\nu}\left|\begin{aligned} \boldsymbol x \newline\boldsymbol y \end{aligned}\right|
-\frac{\hat{s}(\boldsymbol y)}{V_{\Lambda}},
$$

where $V_{\Lambda}=|\det A|$ is the volume of the elementary lattice cell, and the singularity $\hat{s}_\nu(\boldsymbol \cdot)$ is defined as in [epsteinlib](https://github.com/epsteinlib/epsteinlib).

The regularized Epstein zeta function is implemented as
```julia
epsteinzetareg(ν::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})::Complex{Float64}
```
and with optional keyword arguments as
```julia
epsteinzetareg(ν; d, A, x, y)
```
Defaults for `x`, `y`, and `A` are identical to those used in `epsteinzeta`.

## Installation and usage

The library can be installed via

```julia
using Pkg; Pkg.add(url="https://github.com/epsteinlib/EpsteinLib.jl")
```

The following example then computes the Madelung constant to machine precision 

```julia
using EpsteinLib, Printf

ν = 1.0
dim = 3
A = [1.0, 0.0, 0.0,
     0.0, 1.0, 0.0,
     0.0, 0.0, 1.0]
x = [0.0, 0.0, 0.0]
y = [0.5, 0.5, 0.5]

# Calculate Madelung constant
madelung = epsteinzeta(ν, dim, A, x, y)

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