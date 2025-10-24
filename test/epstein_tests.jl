using Test
using LinearAlgebra
using EpsteinZetaFunction
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
        1 1/2;
        0 sqrt(3)/2
    ] # hexagonal lattice matrix
    x = [1/10, 2/10]
    y = [3/10, 4/10]
    @test epsteinzeta(ν; d = d) ≈ -1.9216892211799304 atol = 1e-14
    @test epsteinzeta(ν; A = A) ≈ -1.9999940144822623 atol = 1e-14
    @test epsteinzeta(ν; x = x) ≈ 0.24057039785271267 + 1.039425935548863e-19im atol = 1e-14
    @test epsteinzeta(ν; d = d, A = A, x = x) ≈ 0.1719186692933788 + 8.347631999046855e-22im atol = 1e-14
    @test epsteinzeta(ν; y = y) ≈ -1.2108986338197985 + 2.435700050591742e-19im atol = 1e-14
    @test epsteinzeta(ν; A = A, y = y) ≈ -1.1985660480329243 + 1.1405461005324793e-19im atol = 1e-14
    @test epsteinzeta(ν; x = x, y = y) ≈ 0.8830108146701363 - 0.09354849186479881im atol = 1e-14
    @test epsteinzeta(ν; d = d, A = A, x = x, y = y) ≈ 0.8819439608604308 - 0.10322404491724824im atol = 1e-14
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

@testset "Test single Epstein zeta evaluation" begin

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


@testset "Test single Epstein zeta reg evaluation" begin

    nu = 1/2
    A = [
        1 1/2;
        0 sqrt(3)/2
    ]
    x = [1/10, 2/10]
    y = [3/10, 4/10]

    ref = 0.1225562448097732 + 0.4826367446847953im

    @test epsteinzetareg(nu, A, x, y) ≈ ref atol = 1e-16
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


