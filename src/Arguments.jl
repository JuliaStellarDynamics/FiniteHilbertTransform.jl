#=Arguments.jl

bring in command-line arguments

=#

# To parse the command-line arguments
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
        "--linear"
            help     = "Linear response considered: unstable/damped/neutral"
            arg_type = String
            default  = "damped"
            required = true
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
