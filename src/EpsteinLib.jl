module EpsteinLib

using Epsteinlib_jll, LinearAlgebra

export epsteinzeta, epsteinzetareg

"""
    epsteinzeta(ν::Float64,A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64})
Calls the C function `epsteinZeta` from the shared library.
    double complex epsteinZeta(double nu, unsigned int dim, const double *A, const double *x, const double *y);
Approximates
``Z_{\\nu, A}(x, y) = \\sum_{z \\in A \\mathbb{Z}^d, z \\ne x} \frac{e^{-2\\pi i y \\cdot z}}{|x-z|^\\nu}``
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

function cleanuparguments(d, ν, A, x, y)
    ν = convert(Float64, ν)
    if x==nothing && y==nothing && d==nothing && A==nothing
        throw(ArgumentError("Either d, x, y, or A must be specified"))
    end
    if d==nothing
        if x!=nothing
            d = length(x)
        else
            if y!=nothing
                d = length(y)
            else
                d = size(A, 1)
            end
        end
    end
    d = convert(Int64, d)

    if x == nothing
        x = zeros(d)
    else
        if length(x)==d
            x = convert(Vector{Float64}, x)
        else
            throw(ArgumentError("Incompatible size for x"))
        end
    end
    if y == nothing
        y = zeros(d)
    else
        if length(y)==d
            y = convert(Vector{Float64}, y)
        else
            throw(ArgumentError("Incompatible size for y"))
        end
    end

    if A == nothing
        A = Matrix{Float64}(I, d, d)
    elseif size(A) != (d, d)
        throw(ArgumentError("Incompatible size of A"))
    else
        A = convert(Matrix{Float64}, A)
    end

    return (ν, A, x, y)
end

"""

    epsteinzeta(ν; d, x, y, A)
where d, x, y, A are optional. x and y default to zero of size d, and A to the identity matrix of size d.  

Approximatess
``Z_{\\nu, A}(x, y) = \\sum_{z \\in A \\mathbb{Z}^d, z \\ne x} \frac{e^{-2\\pi i y \\cdot z}}{|x-z|^\\nu}``
if the real part of nu is greater than the system dimension, and the meromorphic continuation otherwise.

"""
function epsteinzeta(
    ν::T0;
    d::Union{Integer,Nothing} = nothing,
    x::Union{Vector{T1},Nothing} = nothing,
    y::Union{Vector{T2},Nothing} = nothing,
    A::Union{Matrix{T3},Nothing} = nothing,
)::Complex{Float64} where {T0<:Real,T1<:Real,T2<:Real,T3<:Real}
    ν, A, x, y = cleanuparguments(d, ν, A, x, y)

    return epsteinzeta(ν, A, x, y)
end

"""
    epsteinzetareg(ν::Float64,A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64})
Calls the C function `epsteinZetaReg` from the shared library.
    double complex epsteinZetaReg(double nu, unsigned int dim, const double *A, const double *x, const double *y);
Calculates a regularization of the Epstein zeta function in the second vector argument.
"""
function epsteinzetareg(
    ν::Float64,
    A::Matrix{Float64},
    x::Vector{Float64},
    y::Vector{Float64},
)::Complex{Float64}
    dim = UInt32(size(A, 1))
    A_flat = vec(permutedims(A))
    return @ccall libepstein.epsteinZetaReg(
        ν::Float64,
        dim::UInt32,
        A_flat::Ref{Float64},
        x::Ref{Float64},
        y::Ref{Float64},
    )::Complex{Float64}
end

"""
    epsteinzetareg(ν; d, x, y, A)
where d, x, y, A are optional. x and y default to zero of size d, and A to the identity matrix of size d.  
Calculates a regularization of the Epstein zeta function in the second vector argument.
"""
function epsteinzetareg(
    ν::T0;
    d::Union{Integer,Nothing} = nothing,
    x::Union{Vector{T1},Nothing} = nothing,
    y::Union{Vector{T2},Nothing} = nothing,
    A::Union{Matrix{T3},Nothing} = nothing,
)::Complex{Float64} where {T0<:Real,T1<:Real,T2<:Real,T3<:Real}
    ν, A, x, y = cleanuparguments(d, ν, A, x, y)

    return epsteinzetareg(ν, A, x, y)
end

end # module
