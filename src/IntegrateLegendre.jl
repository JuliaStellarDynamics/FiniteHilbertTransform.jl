
# To have access to the nodes and weights of the G-L (Gauss-Legendre) quadrature
using FastGaussQuadrature

# Bring in the prefactors
include("Legendre/PrecomputeLegendre.jl")

# Bring in the integration tools
include("Legendre/Legendre.jl")


"""ComputeALegendre

for all values of u:
compute a_k(u) by looping over Legendre weights w(u), P_k(u) values, and G(u) values.

parallel or non-parallel options
"""
function ComputeALegendre(tabwGLquad::Vector{Float64},
                          tabG::Vector{Float64},
                          tabPGLquad::Matrix{Float64},
                          tabINVcGLquad::Vector{Float64})

    K_u = size(tabwGLquad, 1)

    taba = zeros(Float64,K_u)

    Threads.@threads for k=1:(K_u) # Loop over the Legendre functions
        res = 0.0 # Initialisation of the result
        for i=1:K_u # Loop over the G-L nodes
            w = tabwGLquad[i] # Current weight
            G = tabG[i] # Current value of G[u_i]
            P = tabPGLquad[k,i] # Current value of P_k. ATTENTION, to the shift of the array. ATTENTION, to the order of the arguments.
            res += w*G*P # Update of the sum
        end
        res *= tabINVcGLquad[k] # Multiplying by the prefactor. ATTENTION, to the shift of the array

        taba[k] = res # Filling in taba. ATTENTION, to the shift of the array
    end

    return taba

end



"""ComputeIminusXi

wrapper to parallelise calculations of I-Xi

"""
function ComputeIminusXi(tabomega::Vector{Complex{Float64}},
                         taba::Vector{Float64},
                         xmax::Float64,
                         struct_tabLeg::Vector{struct_tabLeg_type})

    # get constants
    K_u    = size(taba,1)
    nomega = size(tabomega,1)

    # define a table to store the value of det[I-Xi].
    tabIminusXi = zeros(Complex{Float64},nomega)

    # loop over all the considered COMPLEX frequencies
    Threads.@threads for iomega=1:nomega

        # ID of the current thread
        thr = Threads.threadid()

        # compute I-Xi(omg) using the parallel containers
        val = GetLegendreIminusXiPlasma(tabomega[iomega],taba,xmax,struct_tabLeg[thr])

        # fill in tabIminusXi
        tabIminusXi[iomega] = val

    end

    return tabIminusXi

end

"""setup_legendre_integration

build various tables for integrating the plasma problem with Legendre

"""
function setup_legendre_integration(K_u::Int64,qself::Float64,xmax::Float64,PARALLEL::Bool=false)

    # Filling in the arrays used in the G-L quadrature (src/Precompute.jl)
    tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = tabGLquad(K_u)

    # compute the function G(u)
    tabG = CGFuncPlasma(tabuGLquad,qself,xmax)

    # compute the coefficients for integration
    taba = ComputeALegendre(tabwGLquad,tabG,tabPGLquad,tabINVcGLquad)

    # set up the table for integration
    struct_tabLeg = initialize_struct_tabLeg(K_u,PARALLEL)

    return taba,struct_tabLeg

end



"""

routine to test the different integration regimes for the plasma case

"""
function test_ninepointsL(taba::Vector{Float64},
                          xmax::Float64,
                          struct_tabLeg::struct_tabLeg_type,
                          digits::Int64=4)
    #=function to test nine unique points for values of D_k

    =#
    upperleft  = -1.5 + 1.5im
    uppercen   =  0.0 + 1.5im
    upperright =  1.5 + 1.5im
    midleft    = -1.5 + 0.0im
    midcen     =  0.0 + 0.0im
    midright   =  1.5 + 0.0im
    lowerleft  = -1.5 - 1.5im
    lowercen   =  0.0 - 1.5im
    lowerright =  1.5 - 1.5im

    println(round(GetLegendreIminusXiPlasma(upperleft,taba,xmax,struct_tabLeg,"unstable"),digits=digits)," || ",
            round(GetLegendreIminusXiPlasma(uppercen,taba,xmax,struct_tabLeg,"unstable"),digits=digits)," || ",
            round(GetLegendreIminusXiPlasma(upperright,taba,xmax,struct_tabLeg,"unstable"),digits=digits))
    println(round(GetLegendreIminusXiPlasma(midleft,taba,xmax,struct_tabLeg,"neutral"),digits=digits)," || ",
            round(GetLegendreIminusXiPlasma(midcen,taba,xmax,struct_tabLeg,"neutral"),digits=digits)," || ",
            round(GetLegendreIminusXiPlasma(midright,taba,xmax,struct_tabLeg,"neutral"),digits=digits))
    println(round(GetLegendreIminusXiPlasma(lowerleft,taba,xmax,struct_tabLeg,"damped"),digits=digits)," || ",
            round(GetLegendreIminusXiPlasma(lowercen,taba,xmax,struct_tabLeg,"damped"),digits=digits)," || ",
            round(GetLegendreIminusXiPlasma(lowerright,taba,xmax,struct_tabLeg,"damped"),digits=digits))
end
