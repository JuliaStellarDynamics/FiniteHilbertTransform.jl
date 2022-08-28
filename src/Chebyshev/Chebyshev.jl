


function get_Chebyshev_Xi(omg::Complex{Float64},
                          taba::Vector{Float64},
                          LINEAR::String="damped";
                          verbose::Int64=0)
    #=
     Defining the correct computation function
    =#

    if (imag(omg) < 0.0)
        if verbose > 2
            println("FiniteHilbertTransform.Chebyshev.get_Chebyshev_Xi: Using DAMPED Chebyshev integration.")
        end

        get_Xi_DAMPED(omg,taba)

    elseif (imag(omg) == 0.0)
        if verbose > 2
            println("FiniteHilbertTransform.Chebyshev.get_Chebyshev_Xi: Using NEUTRAL Chebyshev integration.")
        end

        get_Xi_NEUTRAL(omg,taba)

    else
        # by default use unstable integration
        if verbose > 2
            println("FiniteHilbertTransform.Chebyshev.get_Chebyshev_Xi: Using UNSTABLE Chebyshev integration.")
        end

        get_Xi_UNSTABLE(omg,taba)
    end
end



function get_sumT(omg::Complex{Float64},
                  taba::Vector{Float64})
    #=
    # Computes the sum
    # S_T(\omega) = pi \sum_{k=0}^{K_u-1} a_k T_{k+1}(\omega)
    # using Clenshaw's algorithm
    # ATTENTION, to the range of the sum
    =#
    K_u = size(taba,1)


    b, bp = 0.0+0.0*im, 0.0+0.0*im # Initialising the counters, b_{i} and b_{i+1}
    #####
    # Loop over the terms
    # ATTENTION, this is in decreasing order
    # ATTENTION, the last element is also scrolled over
    for i=K_u:-1:1
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


function get_sumU(omg::Complex{Float64},
                  taba::Vector{Float64})
    #=
    # Computes the sum
    # S_U(\omega) = pi \sum_{k=0}^{K_u-1} a_k U_{k}(\omega)
    # using Clenshaw's algorithm

    =#
    #####

    K_u = size(taba,1)

    b, bp = 0.0+0.0*im, 0.0+0.0*im # Initialising the counters
    #####
    # Loop over the terms
    # ATTENTION, this is in decreasing order
    # ATTENTION, the last element is not scrolled over
    for i=K_u:-1:2
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
