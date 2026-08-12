module EpsteinLib

using Epsteinlib_jll, LinearAlgebra

export epsteinzeta, epsteinzetareg, epsteinzetaaniso, epsteinzetaanisoreg

"""
    checkdimensions(A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})
Verifies that `A` is square and that `x` and `y` match its dimension, and
returns the dimension as the `UInt32` expected by the C interface. Guards the
low-level methods, which pass raw pointers to C without any further checks.
"""
function checkdimensions(A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64})
    d = size(A, 1)
    if size(A, 2) != d
        throw(ArgumentError("A must be square"))
    end
    if length(x) != d
        throw(ArgumentError("Incompatible size for x"))
    end
    if length(y) != d
        throw(ArgumentError("Incompatible size for y"))
    end
    return UInt32(d)
end


"""
    cleanuparguments(d, ν, A, x, y)
Resolves the optional arguments of the keyword methods into the concrete types
expected by the low-level methods. Any of `d`, `A`, `x`, `y` may be `nothing`.

The dimension is taken from `d` if given, otherwise from `x`, `y` or `A`, in
that order. Missing vectors default to zero and a missing `A` to the identity
matrix, both of size `d`. Everything else is converted to `Float64` and checked
for consistent sizes.

Returns `(ν, A, x, y)`. Throws an `ArgumentError` if all of `d`, `A`, `x`, `y`
are `nothing`, or if the given sizes disagree.
"""
function cleanuparguments(d, ν, A, x, y)
    ν = convert(Float64, ν)
    if x === nothing && y === nothing && d === nothing && A === nothing
        throw(ArgumentError("Either d, x, y, or A must be specified"))
    end
    if d === nothing
        if x !== nothing
            d = length(x)
        else
            if y !== nothing
                d = length(y)
            else
                d = size(A, 1)
            end
        end
    end
    d = convert(Int64, d)

    if x === nothing
        x = zeros(d)
    else
        if length(x) == d
            x = convert(Vector{Float64}, x)
        else
            throw(ArgumentError("Incompatible size for x"))
        end
    end
    if y === nothing
        y = zeros(d)
    else
        if length(y) == d
            y = convert(Vector{Float64}, y)
        else
            throw(ArgumentError("Incompatible size for y"))
        end
    end

    if A === nothing
        A = Matrix{Float64}(I, d, d)
    elseif size(A) != (d, d)
        throw(ArgumentError("Incompatible size of A"))
    else
        A = convert(Matrix{Float64}, A)
    end

    return (ν, A, x, y)
end


"""
    cleanupalpha(d, α)
Validates the multi-index ``\\alpha \\in \\mathbb{N}_0^d`` and converts it to the
`Vector{UInt32}` expected by the C interface. Negative entries are rejected
explicitly, so that the caller sees an `ArgumentError` rather than an
`InexactError` from the unsigned conversion.
"""
function cleanupalpha(d, α)
    if length(α) != d
        throw(ArgumentError("Incompatible size of α"))
    end
    if any(a -> a < 0, α)
        throw(ArgumentError("α must have non-negative entries"))
    end
    return convert(Vector{UInt32}, α)
end


"""
    epsteinzeta(ν::Float64,A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64})
Calls the C function `epsteinZeta` from the shared library.
    double complex epsteinZeta(double nu, unsigned int dim, const double *A, const double *x, const double *y);
Approximates
``Z_{\\nu, A}(x, y) = \\sum_{z \\in A \\mathbb{Z}^d, z \\ne x} \\frac{e^{-2\\pi i y \\cdot z}}{|x-z|^\\nu}``
if the real part of nu is greater than the system dimension, and the meromorphic continuation otherwise.
"""
function epsteinzeta(
    ν::Float64,
    A::Matrix{Float64},
    x::Vector{Float64},
    y::Vector{Float64},
)::Complex{Float64}
    dim = checkdimensions(A, x, y)
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
``Z_{\\nu, A}(x, y) = \\sum_{z \\in A \\mathbb{Z}^d, z \\ne x} \\frac{e^{-2\\pi i y \\cdot z}}{|x-z|^\\nu}``
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
    dim = checkdimensions(A, x, y)
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


"""
    epsteinzetaaniso(ν::Float64,A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64},α::Vector{UInt32})
Calls the C function `epsteinZetaAniso` from the shared library.
    double complex epsteinZetaAniso(double nu, unsigned int dim, const double *A, const double *x, const double *y, const unsigned int *alpha);
Approximates the anisotropic Epstein zeta function
``Z_{\\Lambda, \\nu, \\alpha}(x, y) = \\sum_{z \\in \\Lambda'} e^{-2\\pi i y \\cdot z} \\frac{(z-x)^{\\alpha}}{|z-x|^{\\nu}}``
with ``\\Lambda = A \\mathbb{Z}^d`` and ``z^{\\alpha} = z_1^{\\alpha_1} \\dots z_d^{\\alpha_d}``,
if the real part of nu is greater than ``d + |\\alpha|``, and the meromorphic
continuation otherwise. Recovers the Epstein zeta function for ``\\alpha = 0``.
"""
function epsteinzetaaniso(
    ν::Float64,
    A::Matrix{Float64},
    x::Vector{Float64},
    y::Vector{Float64},
    α::Vector{UInt32},
)::Complex{Float64}
    dim = checkdimensions(A, x, y)
    if length(α) != dim
        throw(ArgumentError("Incompatible size of α"))
    end
    A_flat = vec(permutedims(A))
    return @ccall libepstein.epsteinZetaAniso(
        ν::Float64,
        dim::UInt32,
        A_flat::Ref{Float64},
        x::Ref{Float64},
        y::Ref{Float64},
        α::Ref{UInt32},
    )::Complex{Float64}
end

"""
    epsteinzetaaniso(ν, α; d, x, y, A)
where d, x, y, A are optional. x and y default to zero of size d, and A to the
identity matrix of size d. If none of d, x, y, A is given, the dimension is
taken from the length of α.
 
Approximates the anisotropic Epstein zeta function
``Z_{\\Lambda, \\nu, \\alpha}(x, y) = \\sum_{z \\in \\Lambda'} e^{-2\\pi i y \\cdot z} \\frac{(z-x)^{\\alpha}}{|z-x|^{\\nu}}``
if the real part of nu is greater than ``d + |\\alpha|``, and the meromorphic
continuation otherwise.
"""
function epsteinzetaaniso(
    ν::T0,
    α::Vector{T4};
    d::Union{Integer,Nothing} = nothing,
    x::Union{Vector{T1},Nothing} = nothing,
    y::Union{Vector{T2},Nothing} = nothing,
    A::Union{Matrix{T3},Nothing} = nothing,
)::Complex{Float64} where {T0<:Real,T1<:Real,T2<:Real,T3<:Real,T4<:Integer}
    if d === nothing && x === nothing && y === nothing && A === nothing
        d = length(α)
    end
    ν, A, x, y = cleanuparguments(d, ν, A, x, y)
    α = cleanupalpha(size(A, 1), α)

    return epsteinzetaaniso(ν, A, x, y, α)
end

"""
    epsteinzetaanisoreg(ν::Float64,A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64},α::Vector{UInt32})
Calls the C function `epsteinZetaAnisoReg` from the shared library.
    double complex epsteinZetaAnisoReg(double nu, unsigned int dim, const double *A, const double *x, const double *y, const unsigned int *alpha);
Calculates a regularization of the anisotropic Epstein zeta function in the
second vector argument,
``Z^{(\\mathrm{reg})}_{\\Lambda, \\nu, \\alpha}(x, y) = e^{2\\pi i x \\cdot y} Z_{\\Lambda, \\nu, \\alpha}(x, y) - \\frac{\\hat{s}^{(\\alpha)}_{\\nu}(y)}{(-2\\pi i)^{|\\alpha|} V_{\\Lambda}}``
for ``y \\ne 0``, continuously extended to ``y = 0``.
"""
function epsteinzetaanisoreg(
    ν::Float64,
    A::Matrix{Float64},
    x::Vector{Float64},
    y::Vector{Float64},
    α::Vector{UInt32},
)::Complex{Float64}
    dim = checkdimensions(A, x, y)
    if length(α) != dim
        throw(ArgumentError("Incompatible size of α"))
    end
    A_flat = vec(permutedims(A))
    return @ccall libepstein.epsteinZetaAnisoReg(
        ν::Float64,
        dim::UInt32,
        A_flat::Ref{Float64},
        x::Ref{Float64},
        y::Ref{Float64},
        α::Ref{UInt32},
    )::Complex{Float64}
end

"""
    epsteinzetaanisoreg(ν, α; d, x, y, A)
where d, x, y, A are optional. x and y default to zero of size d, and A to the
identity matrix of size d. If none of d, x, y, A is given, the dimension is
taken from the length of α.
 
Calculates a regularization of the anisotropic Epstein zeta function in the
second vector argument, continuously extended to ``y = 0``.
"""
function epsteinzetaanisoreg(
    ν::T0,
    α::Vector{T4};
    d::Union{Integer,Nothing} = nothing,
    x::Union{Vector{T1},Nothing} = nothing,
    y::Union{Vector{T2},Nothing} = nothing,
    A::Union{Matrix{T3},Nothing} = nothing,
)::Complex{Float64} where {T0<:Real,T1<:Real,T2<:Real,T3<:Real,T4<:Integer}
    if d === nothing && x === nothing && y === nothing && A === nothing
        d = length(α)
    end
    ν, A, x, y = cleanuparguments(d, ν, A, x, y)
    α = cleanupalpha(size(A, 1), α)

    return epsteinzetaanisoreg(ν, A, x, y, α)
end

end # module
