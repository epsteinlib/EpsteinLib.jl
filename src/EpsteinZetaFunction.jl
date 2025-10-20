module EpsteinZetaFunction

using Epsteinlib_jll, LinearAlgebra

export epsteinzeta

"""
Calls the C function `epsteinZeta` from the shared library.
    double complex epsteinZeta(double nu, unsigned int dim, const double *A, const double *x, const double *y);
Approximatess
``Z_{\\nu; A}(x, y) = \\sum_{z \\in A \\mathbb{Z}^d, z \\ne x} \frac{-e^{2\\pi i y \\cdot z}}{|x-z|^\\nu}``
if the real part of nu is greater than the system dimension, and the meromorphic continuation otherwise.
"""
function epsteinzeta(
    ν::Float64,
    A::Matrix{Float64},
    x::Vector{Float64},
    y::Vector{Float64},
)::Complex{Float64}
    dim = Base.UInt32(size(A, 1))
    return @ccall libepstein.epsteinZeta(
        ν::Float64,
        dim::UInt32,
        A::Ref{Float64},
        x::Ref{Float64},
        y::Ref{Float64},
    )::Complex{Float64}
end

"""
Approximates
``Z_{\nu, Id}(x, 0) = \\sum_{z \\in \\mathbb Z^d, z \\ne x} |x-z|^{-\\nu}``
if the real part of nu is greater than the system dimension, and the meromorphic continuation otherwise.
"""
function epsteinzeta(ν::Float64, x::Vector{Float64})::Complex{Float64}
    A = Matrix{Float64}(I, length(x), length(x))
    y = zeros(Float64, length(x))
    return epsteinzeta(ν, A, x, y)
end

"""
Approximates
``Z_{\nu, Id}(0, 0) = \\sum_{z \\in \\mathbb Z^d, z \\ne 0} |z|^{-\\nu}``
if the real part of nu is greater than the system dimension, and the meromorphic continuation otherwise.
"""
function epsteinzeta(ν::Float64, d::Int64)::Complex{Float64}
    x = zeros(Float64, d)
    A = Matrix{Float64}(I, d, d)
    y = zeros(Float64, d)
    return epsteinzeta(ν, A, x, y)
end

end # module
