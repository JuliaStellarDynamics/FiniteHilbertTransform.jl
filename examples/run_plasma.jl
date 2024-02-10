#=
An example call:

julia --threads 4 run_plasma.jl --Cmode legendre --parallel 1 --K_u 205 --nOmega 801 --nEta 300 --xmax 20 --Omegamin -4.0 --Omegamax 4.0 --Etamin 0.0 --Etamax 3.0

=#

########################################
# Installing the necessary libraries
########################################
using Pkg
Pkg.add("Plots")
Pkg.add("ArgParse")

using FiniteHilbertTransform


include("PlasmaModel.jl")

using Plots
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
            default  = 20.0
        "--verbose"
            help     = "Set the report flag (larger gives more report)"
            arg_type = Int64
            default  = 1
        "--nOmega"
            help     = "Number of real points to compute"
            arg_type = Int64
            default  = 800
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
            default  = 400
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

"""
    print_arguments(parsed_args)

Prints the parsed arguments in key-value format.

# Arguments
- `parsed_args::Dict`: A dictionary containing the parsed arguments.

# Description
`print_arguments` prints each argument and its corresponding value in the `parsed_args` dictionary.
"""
function print_arguments(parsed_args)
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
end




"""
    get_tabomega(tabOmega::Vector{Float64}, tabEta::Vector{Float64})

Constructs the table of omega values (complex frequency) from the specified real and imaginary components.

# Arguments
- `tabOmega::Vector{Float64}`: Vector containing the real components of frequency values.
- `tabEta::Vector{Float64}`: Vector containing the imaginary components of frequency values.

# Returns
- `tabomega::Vector{Complex{Float64}}`: Vector of complex frequency values.

"""
function get_tabomega(tabOmega::Vector{Float64},tabEta::Vector{Float64})

    nOmega = size(tabOmega,1)
    nEta   = size(tabEta,1)
    nomega = nOmega*nEta

    tabomega = zeros(Complex{Float64},nomega)
    
    icount = 1 

    for iOmega=1:nOmega # Loop over the real part of the frequency
        for iEta=1:nEta # Loop over the complex part of the frequency
            tabomega[icount] = tabOmega[iOmega] + im*tabEta[iEta] # Fill the current value of the complex frequency
            icount += 1 # Update the counter
        end
    end

    return tabomega
end





function main()

    # get the parsed arguments
    parsed_args = parse_commandline()

    if parsed_args["verbose"] > 0
        print_arguments(parsed_args)
    end

    PARALLEL = parsed_args["parallel"]
    qself    = parsed_args["qSELF"]
    xmax     = parsed_args["xmax"]
    K_u      = parsed_args["K_u"]


    if (PARALLEL)
        nb_threads = Threads.nthreads() # Total number of threads for parallel runs
        println("Using $nb_threads threads.")
    end

    # set up the array of frequencies
    tabOmega = collect(range(parsed_args["Omegamin"],parsed_args["Omegamax"],length=parsed_args["nOmega"]))
    tabEta = collect(range(parsed_args["Etamin"],parsed_args["Etamax"],length=parsed_args["nEta"]))
    nomega = parsed_args["nOmega"]*parsed_args["nEta"] # Total number of complex frequencies for which the dispersion function is computed.

    # (flat) array of omega values to check
    tabomega = get_tabomega(tabOmega,tabEta)

    if parsed_args["Cmode"] == "chebyshev"
        taba = setup_chebyshev_integration(K_u,qself,xmax,PARALLEL)
        #test_ninepoints(taba)
        #println(taba)
        @time tabIminusXi = ComputeIminusXi(tabomega,taba,xmax)
    end

    if parsed_args["Cmode"] == "legendre"
        taba,struct_tabLeg = setup_legendre_integration(K_u,qself,xmax,PARALLEL)
        #test_ninepoints()
        #println(taba)
        @time tabIminusXi = ComputeIminusXi(tabomega,taba,xmax,struct_tabLeg)
    end


    # Prefix of the directory where the files are dumped
    prefixnamefile = "data/"

    # Name of the file where the data is dumped
    namefile = prefixnamefile*"data_"*parsed_args["Cmode"]*"_Plasma_Ku_"*string(K_u)*
               "_qSELF_"*string(qself)*"_xmax_"*string(xmax)*".hf5"

    # you can save the data by uncommenting this:
    #print("Dumping the data | ")
    #@time FiniteHilbertTransform.dump_tabIminusXi(namefile,tabomega,tabIminusXi) # Dumping the values of det[I-Xi]

    epsilon_real = reshape(real.(tabIminusXi),parsed_args["nEta"],parsed_args["nOmega"])
    epsilon_imag = reshape(imag.(tabIminusXi),parsed_args["nEta"],parsed_args["nOmega"])
    epsilon = abs.(epsilon_real .+ im * epsilon_imag)


    # Plot
    contour(tabOmega,tabEta,log10.(epsilon), levels=10, color=:black, #levels=[-2.0, -1.5, -1.0, -0.5, -0.25, 0.0], 
            xlabel="Re[ω]", ylabel="Im[ω]", xlims=(-4, 4), ylims=(-3, 0),
            clims=(-2, 0), aspect_ratio=:equal, legend=false)
    savefig("plasmademo.png")


end

main()
