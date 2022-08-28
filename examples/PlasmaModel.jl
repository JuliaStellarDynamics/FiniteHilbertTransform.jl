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
