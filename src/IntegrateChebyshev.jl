using FFTW # To perform discrete sine transform


# Bring in the prefactors
include("Chebyshev/PrecomputeChebyshev.jl")

# Bring in the integration tools
include("Chebyshev/Chebyshev.jl")


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



function test_ninepointsC(taba::Vector{Float64},xmax::Float64,digits::Int64=4)
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

    println(round(get_Chebyshev_IminusXi(upperleft,taba,xmax,"unstable"),digits=digits)," || ",
            round(get_Chebyshev_IminusXi(uppercen,taba,xmax,"unstable"),digits=digits)," || ",
            round(get_Chebyshev_IminusXi(upperright,taba,xmax,"unstable"),digits=digits))
    println(round(get_Chebyshev_IminusXi(midleft,taba,xmax,"neutral"),digits=digits)," || ",
            round(get_Chebyshev_IminusXi(midcen,taba,xmax,"neutral"),digits=digits)," || ",
            round(get_Chebyshev_IminusXi(midright,taba,xmax,"neutral"),digits=digits))
    println(round(get_Chebyshev_IminusXi(lowerleft,taba,xmax,"damped"),digits=digits)," || ",
            round(get_Chebyshev_IminusXi(lowercen,taba,xmax,"damped"),digits=digits)," || ",
            round(get_Chebyshev_IminusXi(lowerright,taba,xmax,"damped"),digits=digits))
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
