




"""
    ChebyshevFHT

"""
struct ChebyshevFHT <: AbstractFHT

    name::String         # FHT name (default Chebyshev)
    Ku::Int64           # number of sample points

    tabu::Array{Float64,1}     # u values (sampling points)
    tabw::Array{Float64,1}     # w values (weights at sampling points)
    tabP::Matrix{Float64}     # P_k(u) values (Ku x Ku)
    tabc::Vector{Float64}     # prefactor at each sampling point

    # arrays for the continuation
    tabPLeg::Array{ComplexF64,1} # Static container for tabPLeg
    tabQLeg::Array{ComplexF64,1} # Static container for tabQLeg
    tabDLeg::Array{ComplexF64,1} # Static container for tabDLeg

end


"""
    ChebyshevFHT(Ku[,name])

Create a ChebyshevFHT structure

"""
function ChebyshevFHT(Ku::Int64;name::String="Chebyshev")

    tabu,tabw,tabc,tabP = tabCquad(Ku)

    return ChebyshevFHT(name,Ku,tabu,tabw,tabP,tabc,
                                  zeros(ComplexF64,Ku),zeros(ComplexF64,Ku),zeros(ComplexF64,Ku))

end




function GettabD!(omg::ComplexF64,
                  FHT::ChebyshevFHT;
                  verbose::Int64=0)
    #=
     Defining the correct computation function
    =#

    if (imag(omg) < 0.0)
        if verbose > 2
            println("FiniteHilbertTransform.Chebyshev.get_Chebyshev_Xi: Using DAMPED Chebyshev integration.")
        end

        get_Xi_DAMPED(omg,FHT.taba)

    elseif (imag(omg) == 0.0)
        if verbose > 2
            println("FiniteHilbertTransform.Chebyshev.get_Chebyshev_Xi: Using NEUTRAL Chebyshev integration.")
        end

        get_Xi_NEUTRAL(omg,FHT.taba)

    else
        # by default use unstable integration
        if verbose > 2
            println("FiniteHilbertTransform.Chebyshev.get_Chebyshev_Xi: Using UNSTABLE Chebyshev integration.")
        end

        get_Xi_UNSTABLE(omg,FHT.taba)
    end
end




function getTu(Ku::Int64,tabuCquad::Vector{Float64})

    tabPCquad = zeros(Float64,Ku,Ku)

    for i=1:Ku # Loop over the nodes

        u = tabuCquad[i] # Current value of the node

        # initialise T_0(u)
        v0 = 1.0

        # initialise T_1(u)
        v1 = u

        tabPCquad[1,i] = v0 # Filling in the value of T_0(u). ATTENTION, to the shift in the array index.
        tabPCquad[2,i] = v1 # Filling in the value of T_1(u). ATTENTION, to the shift in the array index.

        for k=2:(Ku-1) # Loop over the index 2 <= k

            v = 2u*v1 - v0 # Bonnet's recc. relation to get P_k(u)

            # Filling in the value P_k(u). ATTENTION, to the shift in the array index for k
            tabPCquad[k+1,i] = v

            # Shifting the temporary variables
            v0, v1 = v1, v
        end

    end

    return tabPCquad
end



function get_sumT(omg::ComplexF64,
                  taba::Vector{Float64})
    #=
    # Computes the sum
    # S_T(\omega) = pi \sum_{k=0}^{Ku-1} a_k T_{k+1}(\omega)
    # using Clenshaw's algorithm
    # ATTENTION, to the range of the sum
    =#
    Ku = size(taba,1)


    b, bp = 0.0+0.0*im, 0.0+0.0*im # Initialising the counters, b_{i} and b_{i+1}
    #####
    # Loop over the terms
    # ATTENTION, this is in decreasing order
    # ATTENTION, the last element is also scrolled over
    for i=Ku:-1:1
        b, bp = taba[i] + 2.0*omg*b - bp, b # Performing the update
    end
    #####
    # Computing the sum
    # ATTENTION, there is no a_0 in this sum
    # ATTENTION, we used the fact that T_1(w) = w
    S_T = omg*b - bp # ATTENTION, there is no a_0 in this recurrence.
    #####
    S_T *= pi # Rescaling by the factor pi
    #####
    return S_T # Output
end



function get_sumU(omg::ComplexF64,
                  taba::Vector{Float64})
    #=
    # Computes the sum
    # S_U(\omega) = pi \sum_{k=0}^{Ku-1} a_k U_{k}(\omega)
    # using Clenshaw's algorithm

    =#
    #####

    Ku = size(taba,1)

    b, bp = 0.0+0.0*im, 0.0+0.0*im # Initialising the counters
    #####
    # Loop over the terms
    # ATTENTION, this is in decreasing order
    # ATTENTION, the last element is not scrolled over
    for i=Ku:-1:2
        b, bp = taba[i] + 2.0*omg*b - bp, b # Performing the update
    end
    #####
    # Computing the sum
    # ATTENTION, there IS an a_0 in this sum
    # ATTENTION, we used the fact that U_1(w) = 2.0*w
    S_U = taba[1] + 2.0*omg*b - bp # ATTENTION, we must account for the first term a_0
    #####
    S_U *= pi # Rescaling by the factor pi
    #####
    return S_U # Output
end



# bring in the individual cases
include("Unstable.jl")
include("Neutral.jl")
include("Damped.jl")



function GetaXi!(FHT::ChebyshevFHT,
                 tabGXi::AbstractVector{Float64},
                 res::Vector{Float64},warnflag::Vector{Float64})

     # Perfoming the discrete sine transform
     taba_temp = FFTW.r2r(tabGXi,FFTW.RODFT10,1)

     for i=1:FHT.Ku
         res[i] = taba_temp[i] / FHT.Ku
     end
 end

"""

@IMPROVE: how do we handle bad G values?
"""
function GetaXi(FHT::ChebyshevFHT,
                tabGXi::Array{Float64})


    # start with no warnings
    warnflag = zeros(FHT.Ku)

    res = zeros(Float64,FHT.Ku)

    GetaXi!(FHT,tabGXi,res,warnflag)

    return res,warnflag

end
