"""
 the plasma model specific components
"""

import FiniteHilbertTransform

"""
 Pre-computes the needed values of G(u) for the specified plasma model
"""
function CGFuncPlasma(uNodes::Vector{Float64},
                      qSELF::Float64,
                      xmax::Float64)


    K_u = size(uNodes,1)
    tabG = zeros(Float64,K_u)

    # Loop over the nodes
    Threads.@threads for i=1:K_u

        # Current node position
        u_i = uNodes[i]

        # Compute the value of G[u_i]
        tabG[i] = get_G(u_i,qSELF,xmax)

    end

    return tabG
end



#include("../src/Integrate.jl")

"""
the G(u) function for a plasma
"""
function get_G(u::Float64,qSELF::Float64,xmax::Float64)
    x = u*xmax # Value of the velocity
    return qSELF/(sqrt(pi))*x*exp(-x^(2)) # Returning the value of G(u)
end



"""GetLegendreIminusXiPlasma

perform the loop calculation a_k*D_k, after computing D_k

Function that computes the values of I-Xi(omg) for a given complex frequency

"""
function GetLegendreIminusXiPlasma(omg::Complex{Float64},
                                   taba::Vector{Float64},
                                   xmax::Float64,
                                   FHT::FiniteHilbertTransform.AbstractFHT)


    # Rescale the COMPLEX frequency
    K_u = size(taba,1)

    varpi = omg/xmax

    # compute the Hilbert-transformed Legendre functions
    FiniteHilbertTransform.GettabD!(varpi,FHT)

    xi = 0.0 + 0.0*im # Initialise xi

    # loop over the Legendre functions
    for k=1:(K_u)

        # add the contribution
        xi += taba[k]*FHT.tabDLeg[k]
    end

    # compute 1.0 - xi
    IminusXi = 1.0 - xi

    return IminusXi # Output
end




function compute_aChebyshev(tabG::Vector{Float64},
                            PARALLEL::Bool=false)
    #=
     Function that pre-computes the Chebyshev projection
     using a FFT
     @IMPROVE -- the memory handling is not good
    =#

    K_u = size(tabG,1)

    taba = zeros(Float64,K_u)

    # Perfoming the discrete sine transform
    taba_temp = FFTW.r2r(tabG,FFTW.RODFT10,1)

    # ATTENTION, we rescale by K_u, to have the correct normalisation
    if (PARALLEL)
        Threads.@threads for i=1:K_u
            taba[i] = taba_temp[i]/K_u
        end
    else
        for i=1:K_u
            taba[i] = taba_temp[i]/K_u
        end
    end

    return taba

end


"""get_Chebyshev_IminusXi

compute ``({\\bf I}-{\\bf \\Xi}(\\omega))`` for a single frequency.
"""
function get_Chebyshev_IminusXi(omg::Complex{Float64},
                                taba::Vector{Float64},
                                xmax::Float64)
    #=
     Function that computes the values of I-Xi(omg)
     for a given complex frequency
    =#
    #####
    varpi = omg/xmax # Rescaled COMPLEX frequency
    #####
    Xi = get_Chebyshev_Xi(varpi,taba)
    #####
    IminusXi = 1.0 - Xi # Computing the value of 1.0 - xi
    #####
    return IminusXi # Output
end



"""ComputeIminusXi

compute ``({\\bf I}-{\\bf \\Xi}(\\omega))`` for all considered frequencies.
"""
function ComputeIminusXi(tabomega::Vector{Complex{Float64}},
                             taba::Vector{Float64},
                             xmax::Float64)

    K_u = size(taba,1)
    nomega = size(tabomega,1)
    tabIminusXi = zeros(Complex{Float64},nomega) # Table to store the value of det[I-Xi].

    Threads.@threads for iomega=1:nomega # Loop over all the considered COMPLEX frequencies

        thr = Threads.threadid() # ID of the current thread

        tabIminusXi[iomega] = get_Chebyshev_IminusXi(tabomega[iomega],taba,xmax) #

    end

    return tabIminusXi

end

"""setup_chebyshev_integration

fill in the arrays for Chebyshev quadrature,
compute G(u),
compute the a coefficients for integration.
"""
function setup_chebyshev_integration(K_u::Int64,
                                     qself::Float64,
                                     xmax::Float64,
                                     PARALLEL::Bool=false)

    # Filling in the arrays used in the Chebyshev quadrature (src/Precompute.jl)
    tabuCheby = get_tabuCheby(K_u)

    # compute the function G(u)
    tabG = CGFuncPlasma(tabuCheby,qself,xmax)

    # compute the coefficients for integration
    taba = compute_aChebyshev(tabG,PARALLEL)

    return taba

end


"""
    compute_tabG(uNodes, Gfun)

Compute the G function values on the u nodes.
"""
function compute_tabG(uNodes::Vector{Float64},Gfun::Function)

    K_u = size(uNodes,1)        # Nodes number
    tabG = zeros(Float64,K_u)   # Values array

    Threads.@threads for i=1:K_u # Loop over the nodes

        # Current node position
        u_i = uNodes[i]

        # Compute the value of G[u_i]
        tabG[i] = Gfun(u_i)
    end

    return tabG
end



"""ComputeALegendre

for all values of u:
compute a_k(u) by looping over Legendre weights w(u), P_k(u) values, and G(u) values.

parallel or non-parallel options
"""
function ComputeALegendre(FHT::FiniteHilbertTransform.LegendreFHT,tabG::Vector{Float64})

    taba, warnflag = FiniteHilbertTransform.GetaXi(FHT,tabG)

end



"""ComputeIminusXi

wrapper to parallelise calculations of I-Xi

"""
function ComputeIminusXi(tabomega::Vector{Complex{Float64}},
                         taba::Vector{Float64},
                         xmax::Float64,
                         struct_tabLeg::Vector{FiniteHilbertTransform.LegendreFHT})

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
function setup_legendre_integration(Ku::Int64,qself::Float64,xmax::Float64,PARALLEL::Bool=false)

    # Filling in the arrays used in the G-L quadrature (src/Precompute.jl)
    FHT = FiniteHilbertTransform.LegendreFHT(Ku)

    # compute the function G(u)
    tabG = CGFuncPlasma(FHT.tabu,qself,xmax)

    # compute the coefficients for integration
    taba,warnflag = ComputeALegendre(FHT,tabG)

    # set up the table for integration
    FHTlist = [deepcopy(FHT) for k=1:Threads.nthreads()]

    return taba,FHTlist

end
