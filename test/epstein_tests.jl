using Test
using LinearAlgebra
using EpsteinZetaFunction
using SpecialFunctions

@testset "Matches Riemann zeta at d=1" begin
    d = 1
    for ν = 2.0:0.5:5.0
        ref = 2 * zeta(ν)
        @test epsteinzeta(ν; d = 1) ≈ ref atol = 1e-6
        @test epsteinzeta(ν; x = [0.0]) ≈ ref atol = 1e-6
        @test epsteinzeta(ν; d = 1, x = [0.0]) ≈ ref atol = 1e-6
        @test epsteinzeta(ν; y = [0.0]) ≈ ref atol = 1e-6
        @test epsteinzeta(ν; d = 1, y = [0.0]) ≈ ref atol = 1e-6
        @test epsteinzeta(ν; x = [0.0], y = [0.0]) ≈ ref atol = 1e-6
        @test epsteinzeta(ν; d = 1, x = [0.0], y = [0.0]) ≈ ref atol = 1e-6
    end
end

@testset "Test errors" begin
    ν = 2.0
    @test_throws ArgumentError epsteinzeta(ν)
    @test_throws ArgumentError epsteinzeta(ν; d = 1, x = [0.0, 0.0])
    @test_throws ArgumentError epsteinzeta(ν; d = 1, y = [0.0, 0.0])
    @test_throws ArgumentError epsteinzeta(ν; x = [0.0], y = [0.0, 0.0])
    @test_throws ArgumentError epsteinzeta(ν; d = 1, x = [0.0, 0.0], y = [0.0, 0.0])
end

@testset "Compare single non-diagonal evaluation with reference value obtained from the C implementation" begin

    nu = 1/2
    A = [
        1 1/2;
        0 sqrt(3)/2
    ]
    x = [1/10, 2/10]
    y = [3/10, 4/10]

    ref = 0.8819439608604308 - 0.10322404491724824im

    @test epsteinzeta(nu, A, x, y) ≈ ref atol = 1e-16
end


truncated_power(ν, s) = abs(s) > 1e-10 ? abs(s)^(-ν) : 0.0
sum_vertices(ν, N) = sum(truncated_power(ν, sqrt(x^2 + y^2)) for x = (-N):N, y = (-N):N)

@testset "Matches sum at d=2" begin
    d = 2
    N = 5_000
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
    @test z1 == expected

    d = 3
    z2 = epsteinzeta(ν; d = d)
    @test isa(z2, Complex{Float64})

    A2 = Matrix{Float64}(I, d, d)
    x2 = zeros(Float64, d)
    y2 = zeros(Float64, d)
    expected2 = epsteinzeta(ν, A2, x2, y2)
    @test z2 == expected2

end


