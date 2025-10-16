module EpsteinZetaFunction

using LinearAlgebra

export epsteinZeta

# Locate and load the library product

if occursin("arm64-apple-darwin", Sys.MACHINE)
    
    using Libdl
    libepstein_path = joinpath(@__DIR__, "lib", "libepstein.dylib")
    """
    Calls the C function `epsteinZeta` from the shared library.
        double complex epsteinZeta(double nu, unsigned int dim, const double *A, const double *x, const double *y);
    Approximates
    ``Z_{\\nu; A}(x, y) = \\sum_{z \\in A \\mathbb{Z}^d, z \\ne x} \frac{e^{2\\pi i y \\cdot z}}{|x-z|^\\nu}``
    """
    function epsteinZeta(nu::Float64, A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})::Complex{Float64}
        dim = Base.UInt32(size(A, 1))
        return @ccall libepstein_path.epsteinZeta(nu::Float64, dim::UInt32, A::Ref{Float64}, x::Ref{Float64}, y::Ptr{Float64})::Complex{Float64}
    end
else
    error("Unsupported architecture: $(Sys.MACHINE). This library currently supports only arm64-apple-darwin. Waiting for the jll package to appear. Because of a linker error the jll does not support arm64-apple-darwin.")
end

"""
Approximates
``Z_{\nu, Id}(x, 0) = \\sum_{z \\in \\mathbb Z^d, z \\ne x} |x-z|^{-\\nu}``
"""
function epsteinZeta(nu::Float64, x::Vector{Float64})::Complex{Float64}
    A = Matrix{Float64}(I, length(x), length(x))
    y = zeros(Float64, length(x))
    return epsteinZeta(nu, A, x, y)
end

"""
Approximates
``Z_{\nu, Id}(0, 0) = \\sum_{z \\in \\mathbb Z^d, z \\ne 0} |z|^{-\\nu}``
"""
function epsteinZeta(nu::Float64, d::Int64)::Complex{Float64}
    x = zeros(Float64, d)
    A = Matrix{Float64}(I, d, d)
    y = zeros(Float64, d)
    return epsteinZeta(nu, A, x, y)
end

end # module
