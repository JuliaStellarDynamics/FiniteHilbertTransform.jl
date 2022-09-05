"""
 the plasma model specific components
"""

using ArgParse

function parse_commandline()
    #=parse_commandline

    parse command line arguments.
    use --help to see the defaults
    =#
    tabargs = ArgParseSettings()

    @add_arg_table tabargs begin
        "--parallel"
            help     = "Parallel computation: true/false"
            arg_type = Bool
            default  = true
        "--K_u"
            help     = "Number of nodes in the Gauss-Legendre quadrature"
            arg_type = Int64
            default  = 200
        "--qSELF"
            help     = "Self-gravity strength: q < 1 for stable"
            arg_type = Float64
            default  = 0.5
        "--xmax"
            help     = "Truncation range of the velocity range"
            arg_type = Float64
            default  = 5.0
        "--verbose"
            help     = "Set the report flag (larger gives more report)"
            arg_type = Int64
            default  = 1
        "--nOmega"
            help     = "Number of real points to compute"
            arg_type = Int64
            default  = 801
        "--Omegamin"
            help     = "Minimum real frequency"
            arg_type = Float64
            default  = -4.0
        "--Omegamax"
            help     = "Maximum real frequency"
            arg_type = Float64
            default  = 4.0
        "--nEta"
            help     = "Number of imaginary points to compute"
            arg_type = Int64
            default  = 300
        "--Etamin"
            help     = "Minimum imaginary frequency"
            arg_type = Float64
            default  = -3.0
        "--Etamax"
            help     = "Maximum imaginary frequency"
            arg_type = Float64
            default  = -0.01
        "--Cmode"
            help     = "Continuation mode for damped calculation (legendre/chebyshev,rational)"
            arg_type = String
            default  = "legendre"
    end

    return parse_args(tabargs)
end

function print_arguments(parsed_args)
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
end



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



include("../src/Integrate.jl")
include("../src/IntegrateChebyshev.jl")
include("../src/IntegrateLegendre.jl")

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
                                   struct_tabLeg::struct_tabLeg_type)


    # Rescale the COMPLEX frequency
    K_u = size(taba,1)

    varpi = omg/xmax

    # compute the Hilbert-transformed Legendre functions
    get_tabLeg!(varpi,K_u,struct_tabLeg)

    xi = 0.0 + 0.0*im # Initialise xi

    # loop over the Legendre functions
    for k=1:(K_u)

        # add the contribution
        xi += taba[k]*struct_tabLeg.tabDLeg[k]
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
