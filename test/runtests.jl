



using FiniteHilbertTransform
using Test

# parameters common to all tests
Ku = 10
tol = 1e-10

tabu,tabw,tabc,tabP = FiniteHilbertTransform.tabGLquad(Ku)
FHT = FiniteHilbertTransform.LegendreFHT(Ku)


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

@testset "Legendre Unstable Frequency" begin
    @test (@elapsed FiniteHilbertTransform.GettabD!(ϖ,FHT,verbose=4)) < 5e-4
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

@testset "Legendre Neutral Frequency" begin
    @test (@elapsed FiniteHilbertTransform.GettabD!(ϖ,FHT,verbose=4)) < 5e-4
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

@testset "Legendre Damped Frequency" begin
    @test (@elapsed FiniteHilbertTransform.GettabD!(ϖ,FHT,verbose=4)) < 5e-4
    @test (@allocated FiniteHilbertTransform.GettabD!(ϖ,FHT)) == 0
    @test real(FHT.tabDLeg[1]) ≈ -0.0399893282162609 atol=tol
    @test imag(FHT.tabDLeg[1]) ≈ 3.1816033151335357 atol=tol
    @test real(FHT.tabQLeg[1]) == real(FHT.tabDLeg[1])
    @test imag(FHT.tabQLeg[1]) ≈ -3.1015819920460506 atol=tol
    @test real(FHT.tabPLeg[1]) == 1.0
    @test imag(FHT.tabPLeg[1]) == 0.0
end

@testset "Legendre Damped Frequency: Heaviside checks" begin
    ϖ = 0.99 - 0.02im
    FiniteHilbertTransform.GettabD!(ϖ,FHT)
    @test real(FHT.tabDLeg[1]) ≈ -4.48863636973214 atol=tol
    ϖ = 1.0 - 0.02im
    FiniteHilbertTransform.GettabD!(ϖ,FHT)
    @test real(FHT.tabDLeg[1]) ≈ -4.605220183488258 atol=tol
    ϖ = 1.01 - 0.02im
    FiniteHilbertTransform.GettabD!(ϖ,FHT)
    @test real(FHT.tabDLeg[1]) ≈ -4.498635453116724 atol=tol
    ϖ = -1.01 - 0.02im
    FiniteHilbertTransform.GettabD!(ϖ,FHT)
    @test real(FHT.tabDLeg[1]) ≈ 4.498635453116724 atol=tol
end

# check the quadrature approximation
# by checking with all ones, we know the analytic value
@testset "Legendre Quadrature" begin
    # the quadrature values
    Gvals = ones(FHT.Ku)
    @test FiniteHilbertTransform.GetaXi(FHT,Gvals)[1][1] ≈ 1.0 atol=tol
    @test FiniteHilbertTransform.GetaXi(FHT,Gvals)[1][2] ≈ 0.0 atol=tol
    # the flag for a warning
    @test FiniteHilbertTransform.GetaXi(FHT,Gvals)[2] == 0
    # test compatibility version directly
    res = zeros(Float64,Ku)
    warnflag = zeros(Float64,Ku)
    @test FiniteHilbertTransform.GetaXi!(FHT,Gvals,res,warnflag)[2][1] == 0
    #
    ϖ = 0.02 + 0.02im
    @test real(FiniteHilbertTransform.GetIminusXi(ϖ,FiniteHilbertTransform.GetaXi(FHT,Gvals)[1],FHT)) ≈ 1.0399893282162609 atol=tol 
    ϖ = 0.02 + 0.00im
    @test real(FiniteHilbertTransform.GetIminusXi(ϖ,FiniteHilbertTransform.GetaXi(FHT,Gvals)[1],FHT)) ≈ 1.0400053346136993 atol=tol 
    ϖ = 0.02 - 0.02im
    @test real(FiniteHilbertTransform.GetIminusXi(ϖ,FiniteHilbertTransform.GetaXi(FHT,Gvals)[1],FHT)) ≈ 1.0399893282162609 atol=tol 
    #
    # now add a NaN
    Gvals[1] = NaN
    @test FiniteHilbertTransform.GetaXi(FHT,Gvals)[1][1] < 1.0
    # the flag for a warning
    #@test FiniteHilbertTransform.GetaXi(FHT,Gvals)[2] == 1
end


# begin Chebyshev tests
FHT = FiniteHilbertTransform.ChebyshevFHT(Ku)
tabu,tabw,tabc,tabP = FiniteHilbertTransform.tabCquad(Ku)

@testset "Chebyshev Initialisation" begin
    @test (@allocated tabu,tabw,tabc,tabP = FiniteHilbertTransform.tabCquad(Ku)) < 5000
    @test (@allocated FHT = FiniteHilbertTransform.ChebyshevFHT(Ku)) < 5000
    @test FHT.Ku == 10
    @test FHT.tabu[1] ≈ 0.9876883405951378 atol=tol
    @test FHT.tabw[1] == 1.0
end


@testset "Chebyshev Unstable Frequency" begin
    ϖ = 0.02 + 0.02im
    Gvals = ones(FHT.Ku)
    FiniteHilbertTransform.GetaXi(FHT,Gvals)
    @test FiniteHilbertTransform.GetaXi(FHT,Gvals)[1][1] ≈ 1.2784906442999322 atol=tol
    @test FiniteHilbertTransform.GetaXi(FHT,Gvals)[1][2] == 0.0
    @test real(FiniteHilbertTransform.get_sumT(ϖ,FHT.taba)) ≈ 0.10271294847884978 atol=tol
    @test real(FiniteHilbertTransform.get_sumU(ϖ,FHT.taba)) ≈ 3.4519844421597705 atol = tol
    @test real(FiniteHilbertTransform.GettabD!(ϖ,FHT)) ≈ -0.09041331659821551 atol = tol
    ϖ = 0.02 + 0.02im
    @test real(FiniteHilbertTransform.GetIminusXi(ϖ,FiniteHilbertTransform.GetaXi(FHT,Gvals)[1],FHT)) ≈ 1.0904133165982155 atol=tol 
    ϖ = 0.02 + 0.00im
    @test real(FiniteHilbertTransform.GetIminusXi(ϖ,FiniteHilbertTransform.GetaXi(FHT,Gvals)[1],FHT)) ≈ 1.101538304553634 atol=tol 
    ϖ = 0.02 - 0.02im
    @test real(FiniteHilbertTransform.GetIminusXi(ϖ,FiniteHilbertTransform.GetaXi(FHT,Gvals)[1],FHT)) ≈ 1.1150125803594841 atol=tol 
end
