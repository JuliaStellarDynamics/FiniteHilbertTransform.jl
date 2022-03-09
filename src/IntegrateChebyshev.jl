using FFTW # To perform discrete sine transform

# Bring in the prefactors
include("Chebyshev/Precompute.jl")

# Bring in the integration tools
include("Chebyshev/Chebyshev.jl")

# Bring in the velocity
include("PlasmaModel.jl")



function compute_tabG(tabuCheby::Vector{Float64},
                      qSELF::Float64,
                      xmax::Float64,
                      PARALLEL::Bool)
    #=
     Pre-computes the needed values of G(u)
    =#

    K_u = size(tabuCheby,1)
    tabG = zeros(Float64,K_u)

    # Loop over the Chebyshev nodes
    if (PARALLEL)

        Threads.@threads for i=1:K_u

            # Current Chebyshev node
            u_i = tabuCheby[i]

            # Computing the value of G[u_i]
            tabG[i] = get_G(u_i,qSELF,xmax)
        end
    else

        for i=1:K_u

            # Current Chebyshev node
            u_i = tabuCheby[i]

             # Computing the value of G[u_i]
            tabG[i] = get_G(u_i,qSELF,xmax)
        end
    end

    return tabG
end


function compute_taba(tabG::Vector{Float64},
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



function get_IminusXi(omg::Complex{Float64},
                      taba::Vector{Float64},
                      xmax::Float64,
                      LINEAR::String)
    #=
    # Function that computes the values of I-Xi(omg)
    # for a given complex frequency
    =#
    #####
    varpi = omg/xmax # Rescaled COMPLEX frequency
    #####
    Xi = get_Xi(varpi,taba,LINEAR)
    #####
    IminusXi = 1.0 - Xi # Computing the value of 1.0 - xi
    #####
    return IminusXi # Output
end


function compute_tabIminusXi(tabomega::Vector{Complex{Float64}},
                             taba::Vector{Float64},
                             xmax::Float64,
                             LINEAR::String,
                             PARALLEL::Bool=true)
    #=
     Function that computes I-Xi(omg)
     for all the considered frequencies
    =#
    K_u = size(taba,1)
    nomega = size(tabomega,1)
    tabIminusXi = zeros(Complex{Float64},nomega) # Table to store the value of det[I-Xi].

    if (PARALLEL)
        Threads.@threads for iomega=1:nomega # Loop over all the considered COMPLEX frequencies

            thr = Threads.threadid() # ID of the current thread

            tabIminusXi[iomega] = get_IminusXi(tabomega[iomega],taba,xmax,LINEAR) #

        end
    else
        for iomega=1:nomega
            tabIminusXi[iomega] = get_IminusXi(tabomega[iomega],taba,xmax,LINEAR)
        end
    end

    return tabIminusXi

end


function setup_chebyshev_integration(K_u::Int64,
                                     qself::Float64,
                                     xmax::Float64,
                                     PARALLEL::Bool=false)

    # Filling in the arrays used in the G-L quadrature (src/Precompute.jl)
    tabuCheby = get_tabuCheby(K_u)

    # compute the function G(u)
    tabG = compute_tabG(tabuCheby,qself,xmax,PARALLEL)

    # compute the coefficients for integration
    taba = compute_taba(tabG,PARALLEL)

    return taba

end
