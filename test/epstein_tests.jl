using Test
using LinearAlgebra
using EpsteinZetaFunction
using SpecialFunctions

@testset "Matches Riemann zeta at d=1" begin
    d = 1
    for ν = 2.0:0.5:5.0
        @test epsteinZeta(ν, d) ≈ 2 * zeta(ν) atol = 1e-6
    end
end

truncated_power(ν, s) = abs(s) > 1e-10 ? abs(s)^(-ν) : 0.0
sum_vertices(ν, N) = sum(truncated_power(ν, sqrt(x^2 + y^2)) for x = -N:N, y = -N:N)

@testset "Matches sum at d=2" begin
    d = 2
    N = 5_000
    for ν = 3.0:1.0:5.0
        @test epsteinZeta(ν, d) ≈ sum_vertices(ν, N) atol = 1e-2
    end
end

@testset "EpsteinZeta convenience wrappers" begin
    ν = 2.5
    x = [0.1, 0.2]

    z1 = epsteinZeta(ν, x)
    @test isa(z1, Complex{Float64})

    A = Matrix{Float64}(I, length(x), length(x))
    y = zeros(Float64, length(x))
    expected = epsteinZeta(ν, A, x, y)
    @test z1 == expected

    d = 3
    z2 = epsteinZeta(ν, d)
    @test isa(z2, Complex{Float64})

    A2 = Matrix{Float64}(I, d, d)
    x2 = zeros(Float64, d)
    y2 = zeros(Float64, d)
    expected2 = epsteinZeta(ν, A2, x2, y2)
    @test z2 == expected2
end
