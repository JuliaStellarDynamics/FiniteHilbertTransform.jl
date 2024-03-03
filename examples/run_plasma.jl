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


"""
Calculate the G(u) function for a plasma: here, a simple Maxwellian.
see Equation E5 of Fouvry & Prunet (2022)

# Arguments
- `u::Float64`: The velocity parameter.
- `qSELF::Float64`: A parameter specific to the plasma.
- `xmax::Float64`: The maximum value of x.

# Returns
- `Float64`: The value of G(u) for the given parameters.

# Example
```julia
get_G(1.0, 0.5, 10.0)
```
"""
function get_G(u::Float64,qSELF::Float64,xmax::Float64)
    x = u*xmax # Value of the velocity
    return qSELF/(sqrt(pi))*x*exp(-x^(2)) 
end



"""
Pre-compute the values of G(u) for the specified plasma model.

# Arguments
- `uNodes::Vector{Float64}`: Vector of nodes for which G(u) needs to be computed.
- `qSELF::Float64`: A parameter specific to the plasma model.
- `xmax::Float64`: The maximum value of x.

# Returns
- `tabG::Vector{Float64}`: Vector containing pre-computed values of G(u) for each node.

# Example
```julia
CGFuncPlasma([0.1, 0.2, 0.3], 0.5, 10.0)
```
"""
function CGFuncPlasma(uNodes::Vector{Float64}, qSELF::Float64, xmax::Float64)
    K_u = size(uNodes, 1)
    tabG = zeros(Float64, K_u)

    # Loop over the nodes
    Threads.@threads for i = 1:K_u
        # Current node position
        u_i = uNodes[i]

        # Compute the value of G[u_i]
        tabG[i] = get_G(u_i, qSELF, xmax)
    end

    return tabG
end



"""
Compute coefficients a_k(u) for all values of u by looping over Legendre weights, Legendre polynomials, and G(u) values.

# Arguments
- `FHT::FiniteHilbertTransform.LegendreFHT`: LegendreFHT structure.
- `tabG::Vector{Float64}`: Vector of G(u) values.

# Returns
- `taba::Vector{Float64}`: Vector of computed coefficients a_k(u).
- `warnflag::Int`: Warning flag.

# Example
```julia
ComputeALegendre(FHT, tabG)
```
"""
function ComputeA(FHT::FiniteHilbertTransform.AbstractFHT,tabG::Vector{Float64})

    taba, warnflag = FiniteHilbertTransform.GetaXi(FHT,tabG)

end



"""
Wrapper to parallelize calculations of I-Xi.

# Arguments
- `tabomega::Vector{Complex{Float64}}`: Vector of complex frequencies.
- `taba::Vector{Float64}`: Vector containing coefficients for integration.
- `xmax::Float64`: The maximum value of x.
- `struct_tabLeg::Vector{FiniteHilbertTransform.LegendreFHT}`: Vector of LegendreFHT structures.

# Returns
- `tabIminusXi::Vector{Complex{Float64}}`: Vector containing the values of det[I-Xi] for each frequency.

# Example
```julia
ComputeIminusXi([1.0+1.0im, 2.0+2.0im], [0.1, 0.2, 0.3], 10.0, [struct_tabLeg_1, struct_tabLeg_2])
```
"""
function ComputeIminusXi(tabomega::Vector{Complex{Float64}},
                         taba::Vector{Float64},
                         xmax::Float64,
                         struct_tabFHT)

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
        val = FiniteHilbertTransform.GetIminusXi(tabomega[iomega]/xmax,taba,struct_tabFHT[thr])

        # fill in tabIminusXi
        tabIminusXi[iomega] = val

    end

    return tabIminusXi

end


"""
Setup various tables for integrating the plasma problem with Legendre.

# Arguments
- `Ku::Int64`: Integer specifying the number of Legendre nodes.
- `qself::Float64`: A parameter specific to the plasma problem.
- `xmax::Float64`: The maximum value of x.
- `PARALLEL::Bool`: Boolean indicating whether to enable parallel computation. Defaults to `false`.

# Returns
- `taba::Vector{Float64}`: Vector containing coefficients for integration.
- `FHTlist::Vector{FiniteHilbertTransform.LegendreFHT}`: Vector containing LegendreFHT structures for integration.

# Example
```julia
setup_legendre_integration(10, 0.5, 100.0, true)
```
"""
function setup_legendre_integration(Ku::Int64, qself::Float64, xmax::Float64, PARALLEL::Bool=false)
    # Filling in the arrays used in the G-L quadrature
    FHT = FiniteHilbertTransform.LegendreFHT(Ku)

    # Compute the function G(u)
    tabG = CGFuncPlasma(FHT.tabu, qself, xmax)

    # Compute the coefficients for integration
    taba, warnflag = ComputeA(FHT, tabG)

    # Set up the table for integration
    FHTlist = [deepcopy(FHT) for k = 1:Threads.nthreads()]

    return taba, FHTlist
end


function setup_chebyshev_integration(Ku::Int64, qself::Float64, xmax::Float64, PARALLEL::Bool=false)
    # Filling in the arrays used in the Chebyshev quadrature
    FHT = FiniteHilbertTransform.ChebyshevFHT(Ku)

    # Compute the function G(u)
    tabG = CGFuncPlasma(FHT.tabu, qself, xmax)

    # Compute the coefficients for integration
    taba, warnflag = ComputeA(FHT, tabG)

    # Set up the table for integration
    FHTlist = [deepcopy(FHT) for k = 1:Threads.nthreads()]

    return taba, FHTlist
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
            default  = 3.0
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
    CMODE    = parsed_args["Cmode"]


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

    if (CMODE=="legendre")
        taba,struct_tabLeg = setup_legendre_integration(K_u,qself,xmax,PARALLEL)
    elseif (CMODE=="chebyshev")
        taba,struct_tabLeg = setup_chebyshev_integration(K_u,qself,xmax,PARALLEL)
    else
        throw(DomainError(CMODE,"Integration type must be legendre or chebyshev."))
    end

    @time tabIminusXi = ComputeIminusXi(tabomega,taba,xmax,struct_tabLeg)

    epsilon_real = reshape(real.(tabIminusXi),parsed_args["nEta"],parsed_args["nOmega"])
    epsilon_imag = reshape(imag.(tabIminusXi),parsed_args["nEta"],parsed_args["nOmega"])
    epsilon = abs.(epsilon_real .+ im * epsilon_imag)


    # Plot
    contour(tabOmega,tabEta,log10.(epsilon), levels=10, color=:black, 
            xlabel="Re[ω]", ylabel="Im[ω]", xlims=(parsed_args["Omegamin"],parsed_args["Omegamax"]), ylims=(parsed_args["Etamin"],parsed_args["Etamax"]),
            clims=(-2, 0), aspect_ratio=:equal, legend=false)
    savefig("plasmademo.png")


end

main()
