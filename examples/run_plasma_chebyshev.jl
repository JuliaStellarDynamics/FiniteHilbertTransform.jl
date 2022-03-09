#=
An example call:

julia --threads 4 run_plasma_chebyshev.jl --linear damped --parallel 1 --K_u 205 --nOmega 801 --nEta 300 --xmax 20

=#


include("../src/Arguments.jl")
include("../src/IntegrateChebyshev.jl")
include("../src/PlasmaModel.jl")
include("../src/IO.jl")

function get_tabomega(tabOmega::Vector{Float64},tabEta::Vector{Float64})

    nOmega = size(tabOmega,1)
    nEta   = size(tabEta,1)
    nomega = nOmega*nEta

    tabomega = zeros(Complex{Float64},nomega)
    icount = 1 # Initialising the counter
    #####
    for iOmega=1:nOmega # Loop over the real part of the frequency
        for iEta=1:nEta # Loop over the complex part of the frequency
            tabomega[icount] = tabOmega[iOmega] + im*tabEta[iEta] # Filling the current value of the complex frequency
            #####
            icount += 1 # Updating the counter
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
    LINEAR   = parsed_args["linear"]
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

    # build
    taba = setup_chebyshev_integration(K_u,qself,xmax,PARALLEL)

    @time tabIminusXi = compute_tabIminusXi(tabomega,taba,xmax,LINEAR)

    prefixnamefile = "../data/" # Prefix of the directory where the files are dumped
    namefile = prefixnamefile*"data_CHEBYSHEV_Plasma_Ku_"*string(K_u)*
               "_qSELF_"*string(qself)*"_xmax_"*string(xmax)*".hf5" # Name of the file where the data is dumped
    ##################################################
    print("Dumping the data | ") # Printing
    @time dump_tabIminusXi(namefile,tabomega,tabIminusXi) # Dumping the values of det[I-Xi]

end

main()
