



using FiniteHilbertTransform
using Test

Ku = 10
tabu,tabw,tabc,tabP = FiniteHilbertTransform.tabGLquad(Ku)
FHT = FiniteHilbertTransform.LegendreFHT(Ku)

tol = 1e-10

@testset "Legendre Initialisation" begin
    @test (@allocated tabu,tabw,tabc,tabP = FiniteHilbertTransform.tabGLquad(Ku)) < 5000
    @test (@allocated FHT = FiniteHilbertTransform.LegendreFHT(Ku)) < 5000
    @test FHT.Ku == 10
    @test FHT.tabu[1] ≈ -0.9739065285171717 atol=tol
    @test FHT.tabw[1] ≈ 0.06667134430868821 atol=tol
end

# select an unstable frequency
ϖ = 0.02 + 0.02im
FiniteHilbertTransform.GettabD!(ϖ,FHT)

@testset "Unstable Frequency" begin
    @test (@elapsed FiniteHilbertTransform.GettabD!(ϖ,FHT)) < 1.e-4
    @test (@allocated FiniteHilbertTransform.GettabD!(ϖ,FHT)) == 0
    @test real(FHT.tabDLeg[1]) ≈ -0.0399893282162609 atol=tol
    @test imag(FHT.tabDLeg[1]) ≈ 3.1015819920460506 atol=tol
    @test real(FHT.tabDLeg[1]) == real(FHT.tabQLeg[1])
    @test imag(FHT.tabDLeg[1]) == imag(FHT.tabQLeg[1])
    @test real(FHT.tabPLeg[1]) == 0.0
    @test imag(FHT.tabPLeg[1]) == 0.0
end

# select an unstable frequency
ϖ = 0.00 + 0.00im
FiniteHilbertTransform.GettabD!(ϖ,FHT)

@testset "Neutral Frequency" begin
    @test (@elapsed FiniteHilbertTransform.GettabD!(ϖ,FHT)) < 1.e-4
    @test (@allocated FiniteHilbertTransform.GettabD!(ϖ,FHT)) == 0
    @test real(FHT.tabDLeg[1]) == 0.0
    @test imag(FHT.tabDLeg[1]) ≈ 3.141592653589793 atol=tol
    @test real(FHT.tabQLeg[1]) == 0.0
    @test imag(FHT.tabQLeg[1]) == 0.0
    @test real(FHT.tabPLeg[1]) == 1.0
    @test imag(FHT.tabPLeg[1]) == 0.0
end

# select a damped frequency
ϖ = 0.02 - 0.02im
FiniteHilbertTransform.GettabD!(ϖ,FHT)

@testset "Damped Frequency" begin
    @test (@elapsed FiniteHilbertTransform.GettabD!(ϖ,FHT)) < 1.e-4
    @test (@allocated FiniteHilbertTransform.GettabD!(ϖ,FHT)) == 0
    @test real(FHT.tabDLeg[1]) ≈ -0.0399893282162609 atol=tol
    @test imag(FHT.tabDLeg[1]) ≈ 3.1816033151335357 atol=tol
    @test real(FHT.tabQLeg[1]) == real(FHT.tabDLeg[1])
    @test imag(FHT.tabQLeg[1]) ≈ -3.1015819920460506 atol=tol
    @test real(FHT.tabPLeg[1]) == 1.0
    @test imag(FHT.tabPLeg[1]) == 0.0
end

