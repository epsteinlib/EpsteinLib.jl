module EpsteinZetaFunction

using Epsteinlib_jll, LinearAlgebra

export epsteinzeta

"""
    epsteinzeta(ν::Float64,A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64})
Calls the C function `epsteinZeta` from the shared library.
    double complex epsteinZeta(double nu, unsigned int dim, const double *A, const double *x, const double *y);
Approximates
``Z_{\\nu, A}(x, y) = \\sum_{z \\in A \\mathbb{Z}^d, z \\ne x} \frac{-e^{2\\pi i y \\cdot z}}{|x-z|^\\nu}``
if the real part of nu is greater than the system dimension, and the meromorphic continuation otherwise.
"""
function epsteinzeta(
    ν::Float64,
    A::Matrix{Float64},
    x::Vector{Float64},
    y::Vector{Float64},
)::Complex{Float64}
    dim = UInt32(size(A, 1))
    A_flat = vec(permutedims(A))
    return @ccall libepstein.epsteinZeta(
        ν::Float64,
        dim::UInt32,
        A_flat::Ref{Float64},
        x::Ref{Float64},
        y::Ref{Float64},
    )::Complex{Float64}
end

"""

    epsteinzeta(ν; d, x, y, A)
where d, x, y, A are optional. x and y default to zero of size d, and A to the identity matrix of size d.  

Approximatess
``Z_{\\nu, A}(x, y) = \\sum_{z \\in A \\mathbb{Z}^d, z \\ne x} \frac{-e^{2\\pi i y \\cdot z}}{|x-z|^\\nu}``
if the real part of nu is greater than the system dimension, and the meromorphic continuation otherwise.

"""
function epsteinzeta(
    ν::Real;
    d::Union{Integer,Nothing} = nothing,
    x::Union{Vector{Real},Nothing} = nothing,
    y::Union{Vector{Real},Nothing} = nothing,
    A::Union{Matrix{Real},Nothing} = nothing,
)::Complex{Float64}
    if x==nothing && y==nothing && d==nothing
        throw(ArgumentError("Either d, x, or y must be specified"))
    end
    if d==nothing
        if x==nothing
            d = length(y)
            x = zeros(d)
        else
            d = length(x)
            if y==nothing
                y = zeros(d)
            elseif length(y) != d
                throw(ArgumentError("x and y must be the same length"))
            end
        end
    else
        if x == nothing
            x = zeros(d)
        elseif length(x) != d
            throw(ArgumentError("x must be of length d"))
        end
        if y == nothing
            y = zeros(d)
        elseif length(y) != d
            throw(ArgumentError("y must be of length d"))
        end
    end
    x = convert(Vector{Float64}, x)
    y = convert(Vector{Float64}, y)

    if A == nothing
        A = Matrix{Float64}(I, d, d)
    elseif size(A) != (d, d)
        throw(ArgumentError("Dimensions mismatch at A"))
    else
        A = convert(Matrix{Float64}, A)
    end

    return epsteinzeta(ν, A, x, y)
end

end # module
