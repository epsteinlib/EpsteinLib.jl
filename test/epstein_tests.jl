using Test
using LinearAlgebra
using SpecialFunctions

@testset "Matches Riemann zeta at d=1" begin
    d = 1
    for ν = 2.0:0.5:5.0
        ref = 2 * zeta(ν)
        @test epsteinzeta(ν; d = 1) ≈ ref atol = 1e-14
    end
end

@testset "Matches 2D zeta values from C implementation" begin
    ν = 1 / 2
    d = 2
    A = [
        1 1/2
        0 sqrt(3)/2
    ] # hexagonal lattice matrix
    x = [1 / 10, 2 / 10]
    y = [3 / 10, 4 / 10]
    @test epsteinzeta(ν; d = d) ≈ -1.9216892211799304 atol = 1e-14
    @test epsteinzeta(ν; A = A) ≈ -1.9999940144822623 atol = 1e-14
    @test epsteinzeta(ν; x = x) ≈ 0.24057039785271267 + 1.039425935548863e-19im atol = 1e-14
    @test epsteinzeta(ν; d = d, A = A, x = x) ≈ 0.1719186692933788 + 8.347631999046855e-22im atol =
        1e-14
    @test epsteinzeta(ν; y = y) ≈ -1.2108986338197985 + 2.435700050591742e-19im atol = 1e-14
    @test epsteinzeta(ν; A = A, y = y) ≈ -1.1985660480329243 + 1.1405461005324793e-19im atol =
        1e-14
    @test epsteinzeta(ν; x = x, y = y) ≈ 0.8830108146701363 - 0.09354849186479881im atol =
        1e-14
    @test epsteinzeta(ν; d = d, A = A, x = x, y = y) ≈
          0.8819439608604308 - 0.10322404491724824im atol = 1e-14
end

@testset "Test single Epstein zeta evaluation" begin

    nu = 1 / 2
    A = [
        1 1/2
        0 sqrt(3)/2
    ]
    x = [1 / 10, 2 / 10]
    y = [3 / 10, 4 / 10]

    ref = 0.8819439608604308 - 0.10322404491724824im

    @test epsteinzeta(nu, A, x, y) ≈ ref atol = 2e-16
end


@testset "Test single Epstein zeta reg evaluation" begin

    nu = 1 / 2
    A = [
        1 1/2
        0 sqrt(3)/2
    ]
    x = [1 / 10, 2 / 10]
    y = [3 / 10, 4 / 10]

    ref = 0.1225562448097732 + 0.4826367446847953im

    @test epsteinzetareg(nu, A, x, y) ≈ ref atol = 2e-16
end


truncated_power(ν, s) = abs(s) > 1e-10 ? abs(s)^(-ν) : 0.0
sum_vertices(ν, N) = sum(truncated_power(ν, sqrt(x^2 + y^2)) for x = (-N):N, y = (-N):N)

@testset "Matches sum at d=2" begin
    d = 2
    N = 1_000
    for ν = 3.0:1.0:5.0
        @test epsteinzeta(ν; d = d) ≈ sum_vertices(ν, N) atol = 1e-2
    end
end


@testset "epsteinzeta convenience wrappers" begin
    ν = 2.5
    x = [0.1, 0.2]

    z1 = epsteinzeta(ν; x = x)
    @test isa(z1, Complex{Float64})

    A = Matrix{Float64}(I, length(x), length(x))
    y = zeros(Float64, length(x))
    expected = epsteinzeta(ν, A, x, y)
    @test z1 ≈ expected

    d = 3
    z2 = epsteinzeta(ν; d = d)
    @test isa(z2, Complex{Float64})

    A_id = Matrix{Float64}(I, d, d)
    zeros_d = zeros(Float64, d)
    expected2 = epsteinzeta(ν, A_id, zeros_d, zeros_d)
    @test z2 ≈ expected2

    @test epsteinzeta(ν; x = zeros_d) ≈ expected2
    @test epsteinzeta(ν; d = d, x = zeros_d) ≈ expected2
    @test epsteinzeta(ν; y = zeros_d) ≈ expected2
    @test epsteinzeta(ν; d = d, y = zeros_d) ≈ expected2
    @test epsteinzeta(ν; x = zeros_d, y = zeros_d) ≈ expected2
    @test epsteinzeta(ν; d = d, x = zeros_d, y = zeros_d) ≈ expected2

    @test epsteinzeta(ν; A = A_id) ≈ expected2
    @test epsteinzeta(ν; d = d, A = A_id) ≈ expected2
    @test epsteinzeta(ν; x = zeros_d, A = A_id) ≈ expected2
    @test epsteinzeta(ν; d = d, x = zeros_d, A = A_id) ≈ expected2
    @test epsteinzeta(ν; y = zeros_d, A = A_id) ≈ expected2
    @test epsteinzeta(ν; d = d, y = zeros_d) ≈ expected2
    @test epsteinzeta(ν; x = zeros_d, y = zeros_d) ≈ expected2
    @test epsteinzeta(ν; d = d, x = zeros_d, y = zeros_d) ≈ expected2
end

@testset "General Real types work" begin
    ν = 2.0
    ref = epsteinzeta(ν; d = 3)

    @test epsteinzeta(2; d = Int32(3)) ≈ ref
    @test epsteinzeta(2; x = [0, 0, 0]) ≈ ref
    @test epsteinzeta(2; y = [0, 0, 0]) ≈ ref
    @test epsteinzeta(2; A = [1 0 0; 0 1 0; 0 0 1]) ≈ ref
end

@testset "Test errors" begin
    ν = 2.0
    @test_throws ArgumentError epsteinzeta(ν)
    @test_throws ArgumentError epsteinzeta(ν; d = 1, x = [0.0, 0.0])
    @test_throws ArgumentError epsteinzeta(ν; d = 1, y = [0.0, 0.0])
    @test_throws ArgumentError epsteinzeta(ν; x = [0.0], y = [0.0, 0.0])
    @test_throws ArgumentError epsteinzeta(ν; d = 1, x = [0.0, 0.0], y = [0.0, 0.0])

    A = Matrix{Float64}(I, 1, 1)
    @test_throws ArgumentError epsteinzeta(ν; d = 2, A = A)
    @test_throws ArgumentError epsteinzeta(ν; x = [0.0, 0.0], A = A)
    @test_throws ArgumentError epsteinzeta(ν; y = [0.0, 0.0], A = A)
end

@testset "Low-level methods validate dimensions" begin
    ν = 2.0
    A3 = Matrix{Float64}(I, 3, 3)

    # x or y shorter than dim — C would read past the end of the array
    @test_throws ArgumentError epsteinzeta(ν, A3, zeros(2), zeros(3))
    @test_throws ArgumentError epsteinzeta(ν, A3, zeros(3), zeros(2))
    @test_throws ArgumentError epsteinzetareg(ν, A3, zeros(2), zeros(3))
    @test_throws ArgumentError epsteinzetareg(ν, A3, zeros(3), zeros(2))

    # x or y longer than dim — trailing entries silently ignored
    @test_throws ArgumentError epsteinzeta(ν, A3, zeros(4), zeros(3))
    @test_throws ArgumentError epsteinzeta(ν, A3, zeros(3), zeros(4))
    @test_throws ArgumentError epsteinzetareg(ν, A3, zeros(4), zeros(3))

    # non-square A — flattening produces a silently wrong matrix
    Awide = Matrix{Float64}(I, 2, 3)
    Atall = Matrix{Float64}(I, 3, 2)
    @test_throws ArgumentError epsteinzeta(ν, Awide, zeros(2), zeros(2))
    @test_throws ArgumentError epsteinzeta(ν, Atall, zeros(3), zeros(3))
    @test_throws ArgumentError epsteinzetareg(ν, Awide, zeros(2), zeros(2))
    @test_throws ArgumentError epsteinzetareg(ν, Atall, zeros(3), zeros(3))

    # the guard must not reject valid input
    @test epsteinzeta(ν, A3, zeros(3), zeros(3)) isa Complex{Float64}
    @test epsteinzetareg(ν, A3, zeros(3), zeros(3)) isa Complex{Float64}

    # the same guards apply to the anisotropic methods
    α3 = UInt32[1, 0, 2]
    @test_throws ArgumentError epsteinzetaaniso(ν, A3, zeros(2), zeros(3), α3)
    @test_throws ArgumentError epsteinzetaaniso(ν, A3, zeros(3), zeros(4), α3)
    @test_throws ArgumentError epsteinzetaanisoreg(ν, A3, zeros(2), zeros(3), α3)
    @test_throws ArgumentError epsteinzetaanisoreg(ν, A3, zeros(3), zeros(4), α3)
    @test_throws ArgumentError epsteinzetaaniso(ν, Awide, zeros(2), zeros(2), UInt32[1, 0])
    @test_throws ArgumentError epsteinzetaanisoreg(ν, Atall, zeros(3), zeros(3), α3)

    # α itself must match the dimension
    @test_throws ArgumentError epsteinzetaaniso(ν, A3, zeros(3), zeros(3), UInt32[1, 0])
    @test_throws ArgumentError epsteinzetaaniso(
        ν,
        A3,
        zeros(3),
        zeros(3),
        UInt32[1, 0, 2, 1],
    )
    @test_throws ArgumentError epsteinzetaanisoreg(ν, A3, zeros(3), zeros(3), UInt32[1, 0])

    # the guard must not reject valid input
    @test epsteinzetaaniso(ν, A3, zeros(3), zeros(3), α3) isa Complex{Float64}
    @test epsteinzetaanisoreg(ν, A3, zeros(3), zeros(3), α3) isa Complex{Float64}
end


@testset "Matches anisotropic values from reference Mathematica wrapper" begin
    A_hex = [
        1 1/2
        0 sqrt(3)/2
    ] # hexagonal lattice matrix
    A_3d = [
        1 1/3 0
        0 1 1/4
        0 0 1
    ]

    # (ν, A, x, y, α, aniso reference, anisoreg reference, atol)
    cases = [
        (
            1 / 2,
            A_hex,
            [1 / 10, 2 / 10],
            [3 / 10, 4 / 10],
            UInt32[1, 2],
            0.07983310242582402 + 0.1423578960983542im,
            -0.029229875472268955 + 0.0011338524512871971im,
            0.0,
        ),
        (
            5 / 2,
            A_3d,
            [1 / 10, 2 / 10, 3 / 10],
            [2 / 5, 1 / 5, 3 / 5],
            UInt32[1, 0, 2],
            0.010647919006801941 + 0.188941143464776im,
            -0.18923689964399026 - 0.17797688215910856im,
            0.0,
        ),
        (
            1 / 2,
            A_hex,
            [1 / 10, 2 / 10],
            [3 / 10, 4 / 10],
            UInt32[4, 1],
            0.10323147765611702 + 0.16604413648941896im,
            -0.026299295350110257 + 0.0072935367074584145im,
            0.0,
        ),
        # α_1 odd with x_1 = y_1 = 0: both vanish. Mathematica returns roundoff
        # at 1e-19 for the regularized case, so compare against exact zero with
        # an absolute tolerance instead.
        (
            1 / 2,
            A_hex,
            [0.0, 2 / 10],
            [0.0, 4 / 10],
            UInt32[1, 2],
            0.0 + 0.0im,
            0.0 + 0.0im,
            1e-14,
        ),
    ]

    @testset "ν=$ν, α=$(Int.(α))" for (ν, A, x, y, α, ref, refreg, tol) in cases
        @testset "$f" for (f, r) in ((epsteinzetaaniso, ref), (epsteinzetaanisoreg, refreg))
            @test f(ν, A, x, y, α) ≈ r rtol = 1e-14 atol = tol
        end
    end
end

@testset "Anisotropic structure" begin
    d = 3
    A = Matrix{Float64}(I, d, d)
    x = zeros(d)
    y = fill(0.5, d)
    ν = 4.0
    zeroα = zeros(UInt32, d)

    # α = 0 must reproduce the isotropic functions exactly
    @test epsteinzetaaniso(ν, A, x, y, zeroα) == epsteinzeta(ν, A, x, y)
    @test epsteinzetaanisoreg(ν, A, x, y, zeroα) == epsteinzetareg(ν, A, x, y)

    # keyword form agrees with the positional form
    @test epsteinzetaaniso(ν, [0, 1, 0]; x = x, y = y, A = A) ==
          epsteinzetaaniso(ν, A, x, y, UInt32[0, 1, 0])

    # dimension inferred from α alone
    @test epsteinzetaaniso(ν, [0, 1, 0]) isa Complex{Float64}

    # integer input is accepted
    @test epsteinzetaaniso(4, [0, 1, 0]; x = [0, 0, 0], y = [1, 1, 1]) isa Complex{Float64}

    @test_throws ArgumentError epsteinzetaaniso(ν, [0, -1, 0]; A = A)
end
