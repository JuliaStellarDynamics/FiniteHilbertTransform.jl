
# To have access to the nodes and weights of the G-L (Gauss-Legendre) quadrature
using FastGaussQuadrature

# Bring in the prefactors
include("Legendre/Precompute.jl")

# Bring in the integration tools
include("Legendre/Legendre.jl")

include("PlasmaModel.jl")

# everything from here down is Legendre specific right now.
function compute_tabG(tabuGLquad::Vector{Float64},
                      qSELF::Float64,
                      xmax::Float64,
                      PARALLEL::Bool=false)
    #=

    ##################################################
    # Function that pre-computes the values of tabG
    ##################################################
    =#

    # how many u values are we evaluating?

    K_u = size(tabuGLquad, 1)

    tabG = zeros(Float64,K_u)

    if (PARALLEL) # The calculation is made in parallel
        Threads.@threads for i=1:K_u # Loop over the G-L nodes
            u_i = tabuGLquad[i] # Current value of the G-L node
            #####
            tabG[i] = get_G(u_i,qSELF,xmax) # Computing the value of G[u_i]
        end
    else # The calculation is not made in parallel
        for i=1:K_u # Loop over the G-L nodes
            u_i = tabuGLquad[i] # Current value of the G-L node
            #####
            tabG[i] = get_G(u_i,qSELF,xmax) # Computing the value of G[u_i]
        end
    end

    # no return; tabG is modified in place.
    return tabG
end


function compute_taba(tabwGLquad::Vector{Float64},
              tabG::Vector{Float64},
              tabPGLquad::Matrix{Float64},
              tabINVcGLquad::Vector{Float64},
              PARALLEL::Bool)
    #=##################################################
    # Function that pre-computes the coefficients a
    # using the G-L quadrature
    ##################################################
    =#
    K_u = size(tabwGLquad, 1)

    taba = zeros(Float64,K_u)


    if (PARALLEL) # The calculation is made in parallel
        Threads.@threads for k=0:(K_u-1) # Loop over the Legendre functions
            res = 0.0 # Initialisation of the result
            for i=1:K_u # Loop over the G-L nodes
                w = tabwGLquad[i] # Current weight
                G = tabG[i] # Current value of G[u_i]
                P = tabPGLquad[k+1,i] # Current value of P_k. ATTENTION, to the shift of the array. ATTENTION, to the order of the arguments.
                res += w*G*P # Update of the sum
            end
            res *= tabINVcGLquad[k+1] # Multiplying by the prefactor. ATTENTION, to the shift of the array
            #####
            taba[k+1] = res # Filling in taba. ATTENTION, to the shift of the array
        end
    else # The calculation is not made in parallel
        for k=0:(K_u-1) # Loop over the Legendre functions
            res = 0.0 # Initialisation of the result
            for i=1:K_u # Loop over the G-L nodes
                w = tabwGLquad[i] # Current weight
                G = tabG[i] # Current value of G[u_i]
                P = tabPGLquad[k+1,i] # Current value of P_k. ATTENTION, to the shift of the array. ATTENTION, to the order of the arguments.
                res += w*G*P # Update of the sum
            end
            res *= tabINVcGLquad[k+1] # Multiplying by the prefactor. ATTENTION, to the shift of the array
            #####
            taba[k+1] = res # Filling in taba. ATTENTION, to the shift of the array
        end
    end

    return taba

end



function get_IminusXi(omg::Complex{Float64},
                      taba::Vector{Float64},
                      xmax::Float64,
                      K_u::Int64,
                      struct_tabLeg::struct_tabLeg_type,
                      LINEAR::String="damped")
    #=
    ##################################################
    # Function that computes the values of I-Xi(omg)
    # for a given complex frequency
    ##################################################
    =#
    #####
    varpi = omg/(xmax) # Rescaled COMPLEX frequency
    #####
    get_tabLeg!(varpi,K_u,struct_tabLeg,LINEAR) # Computing the Hilbert-transformed Legendre functions
    tabDLeg = struct_tabLeg.tabDLeg # Name of the array where the D_k(w) are stored
    #####
    xi = 0.0 + 0.0*im # Initialisation of the xi
    #####
    for k=0:(K_u-1) # Loop over the Legendre functions
        xi += taba[k+1]*tabDLeg[k+1] # Adding a contribution. ATTENTION, to the shift of the array.
    end
    #####
    IminusXi = 1.0 - xi # Computing the value of 1.0 - xi
    return IminusXi # Output
end




function compute_tabIminusXi(tabomega::Vector{Complex{Float64}},
                             taba::Vector{Float64},
                             xmax::Float64,
                             #K_u::Int64,
                             struct_tabLeg::Vector{struct_tabLeg_type},
                             LINEAR::String="damped")
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

        val = get_IminusXi(tabomega[iomega],taba,xmax,K_u,struct_tabLeg[thr]) #

        #Computing I-Xi(omg) using the parallel containers

        tabIminusXi[iomega] = val # Filling in tabIminusXi
    end

    return tabIminusXi

end


function setup_legendre_integration(K_u::Int64,qself::Float64,xmax::Float64,PARALLEL::Bool=false)

    # Filling in the arrays used in the G-L quadrature (src/Precompute.jl)
    tabuGLquad,tabwGLquad,tabINVcGLquad,tabPGLquad = tabGLquad(K_u)

    # compute the function G(u)
    tabG = compute_tabG(tabuGLquad,qself,xmax,PARALLEL)

    # compute the coefficients for integration
    taba = compute_taba(tabwGLquad,tabG,tabPGLquad,tabINVcGLquad,PARALLEL)

    # set up the table for integration
    struct_tabLeg = initialize_struct_tabLeg(K_u,PARALLEL)

    return taba,struct_tabLeg

end
