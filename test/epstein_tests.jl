using Test
using LinearAlgebra
using EpsteinZetaFunction
using SpecialFunctions

@testset "Matches Riemann zeta at d=1" begin
    d = 1
    for nu = 2.0:0.5:5.0 
        @test epsteinZeta(nu, d) ≈ 2*zeta(nu) atol=1e-6
    end
end

@testset "EpsteinZeta convenience wrappers" begin
    nu = 2.5
    x = [0.1, 0.2]

    z1 = epsteinZeta(nu, x)
    @test isa(z1, Complex{Float64})    

    A = Matrix{Float64}(I, length(x), length(x))
    y = zeros(Float64, length(x))
    expected = epsteinZeta(nu, A, x, y)
    @test z1 == expected

    d = 3
    z2 = epsteinZeta(nu, d)
    @test isa(z2, Complex{Float64})

    A2 = Matrix{Float64}(I, d, d)
    x2 = zeros(Float64, d)
    y2 = zeros(Float64, d)
    expected2 = epsteinZeta(nu, A2, x2, y2)
    @test z2 == expected2
end
