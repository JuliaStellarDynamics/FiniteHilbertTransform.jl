
# To have access to the nodes and weights of the G-L (Gauss-Legendre) quadrature
using FastGaussQuadrature

# Bring in the prefactors
include("Legendre/PrecomputeLegendre.jl")

# Bring in the integration tools
include("Legendre/Legendre.jl")


"""compute_aLegendre

for all values of u:
compute a_k(u) by looping over Legendre weights w(u), P_k(u) values, and G(u) values.

parallel or non-parallel options
"""
function compute_aLegendre(tabwGLquad::Vector{Float64},
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


"""GetLegendreIminusXiPlasma

perform the loop calculation a_k*D_k, after computing D_k

Function that computes the values of I-Xi(omg) for a given complex frequency

"""
function GetLegendreIminusXiPlasma(omg::Complex{Float64},
                                   taba::Vector{Float64},
                                   xmax::Float64,
                                   struct_tabLeg::struct_tabLeg_type,
                                   LINEAR::String="damped")


    # Rescale the COMPLEX frequency
    K_u = size(taba,1)

    varpi = omg/xmax

    # compute the Hilbert-transformed Legendre functions
    get_tabLeg!(varpi,K_u,struct_tabLeg,LINEAR)

    xi = 0.0 + 0.0*im # Initialise xi

    # loop over the Legendre functions
    for k=1:(K_u)

        # add the contribution
        xi += taba[k]*struct_tabLeg.tabDLeg[k]
    end

    IminusXi = 1.0 - xi # Compute 1.0 - xi
    return IminusXi # Output
end


"""compute_tabIminusXi

wrapper to parallelise calculations of I-Xi

"""
function compute_tabIminusXi(tabomega::Vector{Complex{Float64}},
                             taba::Vector{Float64},
                             xmax::Float64,
                             struct_tabLeg::Vector{struct_tabLeg_type},
                             LINEAR::String)
    #=
     Function that computes I-Xi(omg)
     for all the considered frequencies
     ATTENTION, the parallelism is hard-coded
    =#
    K_u = size(taba,1)
    nomega = size(tabomega,1)
    tabIminusXi = zeros(Complex{Float64},nomega) # Table to store the value of det[I-Xi].

    Threads.@threads for iomega=1:nomega # Loop over all the considered COMPLEX frequencies
        #####
        thr = Threads.threadid() # ID of the current thread

        val = GetLegendreIminusXiPlasma(tabomega[iomega],taba,xmax,struct_tabLeg[thr],LINEAR) #

        #Computing I-Xi(omg) using the parallel containers

        tabIminusXi[iomega] = val # Filling in tabIminusXi
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
    tabG = compute_tabG(tabuGLquad,qself,xmax)

    # compute the coefficients for integration
    taba = compute_aLegendre(tabwGLquad,tabG,tabPGLquad,tabINVcGLquad)

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
