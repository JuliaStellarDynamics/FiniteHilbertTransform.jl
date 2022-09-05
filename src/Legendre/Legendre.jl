"""Legendre.jl
 Functions that compute the dispersion function D_k(w)

 These are defined as
 D_k(w) = INT[P_k(u)/(u-w),{u,-1,1}],
 where the integral follows Landau's prescription

 As such, one has
 D_k(w) = Q_k(w)                                   for Im[w] > 0
        = Q_k(w) +   i*pi*P_k(w)*Heaviside[w]      for Im[w] = 0
        = Q_k(w) + 2*i*pi*P_k(w)*Heaviside[Re[w]]  for Im[w] < 0

 Here, the three different cases correspond
 respectively to unstable; neutral; damped modes

"""

"""
    structLegendreFHTtype

"""
struct structLegendreFHTtype

    name::String         # FHT name (default Legendre)
    Ku::Int64           # number of sample points

    tabu::Array{Float64,1}     # u values (sampling points)
    tabw::Array{Float64,1}     # w values (weights at sampling points)
    tabP::Matrix{Float64}     # P_k(u) values (Ku x Ku)
    tabc::Vector{Float64}     # prefactor at each sampling point

    # arrays for the continuation
    tabPLeg::Array{Complex{Float64},1} # Static container for tabPLeg
    tabQLeg::Array{Complex{Float64},1} # Static container for tabQLeg
    tabDLeg::Array{Complex{Float64},1} # Static container for tabDLeg

end


"""
    LegendreFHTcreate(Ku[name, dimension, lmax, nmax, G, rb])

Create a structLegendreFHTtype structure

"""
function LegendreFHTcreate(Ku::Int64;name::String="Legendre")

    tabu,tabw,tabc,tabP = tabGLquad(Ku)

    return structLegendreFHTtype(name,Ku,tabu,tabw,tabP,tabc,zeros(Complex{Float64},Ku),zeros(Complex{Float64},Ku),zeros(Complex{Float64},Ku))

end



"""initialize_struct_tabLeg
parallelise initialisation of struct_tabLeg_type
"""
function initialize_struct_tabLeg(Ku::Int64,PARALLEL::Bool)

    if (PARALLEL)
        # For parallel runs, we also create containers for each thread
        # Basis container for every threads
        nb_threads = Threads.nthreads()
        struct_tabLeg = [struct_tabLeg_create(Ku) for thr=1:nb_threads]
    else
        # Create a default struct_tabLeg used in serial runs
        # Create one container for non-parallel runs
        struct_tabLeg = struct_tabLeg_create(Ku)
    end

    return struct_tabLeg

end


"""
fill struct_tabLeg at a given complex frequency for the integration being considered

integration style selection is automatic: if you want to specify a type, call out to the specific integration method.

"""
function get_tabLeg!(omg::Complex{Float64},
                     struct_tabLeg::structLegendreFHTtype;
                     verbose::Int64=0)

    # check the imaginary sign. if negative, use damped integration
    if (imag(omg) < 0.0)
        if verbose > 2
            println("FiniteHilbertTransform.Legendre.get_tabLeg!: Using DAMPED Legendre integration.")
        end

        tabLeg!_DAMPED(omg,struct_tabLeg)

    # if exactly zero, use neutral mode calculation
    elseif (imag(omg) == 0.0)
        if verbose > 2
            println("FiniteHilbertTransform.Legendre.get_tabLeg!: Using NEUTRAL Legendre integration.")
        end

        tabLeg!_NEUTRAL(omg,struct_tabLeg)

    # by default, use unstable integration
    else

        if verbose > 2
            println("FiniteHilbertTransform.Legendre.get_tabLeg!: Using UNSTABLE Legendre integration.")
        end

        tabLeg!_UNSTABLE(omg,struct_tabLeg)
    end
end


"""Heaviside
Heaviside function on the interval [-1,1]
Here, H[x] has a REAL argument, and returns
0   for x < -1
1/2 for x = -1
1   for -1 < x < 1
1/2 for x = 1
0   for 1 < x
ATTENTION, the equality tests on Float64 might not be very robust
"""
function Heaviside(x::Float64)
    if     (x <  -1.0)      # Left of the interval
        return 0.0
    elseif (x == -1.0)      # Left edge of the interval
        return 0.5
    elseif (-1.0 < x < 1.0) # Within the interval
        return 1.0
    elseif (x ==  1.0)      # Right edge of the interval
        return 0.5
    elseif (1.0 < x)        # Right of the interval
        return 0.0
    end
end



function tabQLeg!(omg::Complex{Float64},
                  val_0::Complex{Float64},
                  val_1::Complex{Float64},
                  tabQLeg::Array{Complex{Float64},1})
    #=tabQLeg
     Function that pre-computes the Hilbert-transformed
     Legendre functions for a given complex frequency
     Q_k(w) = INT[P_k(u)/(u-w),{u,-1,1}]
     ATTENTION, this is == -2 q_k(w)
     with q_k(w) the Legendre functions of the second kind
     for REAL values of w.

     Arguments are:
     + omg: COMPLEX frequency.      ATTENTION, has to be complex.
     + val_0: Initial value in k=0. ATTENTION, has to be complex.
     + val_1: Initial value in k=1. ATTENTION, has to be complex.
     + tabQLeg: Container where to store the results
     There are two possible algorithms:
     + Close to the real line [-1,1], we use an upward recurrence
     + Far-away from the real line[-1,1], we use a downward recurrence
     The transition from the two regimes follows from the thesis
     Stable Implementation of Three-Term Recurrence Relations
     Pascal Frederik Heiter, June, 2010
     https://www.uni-ulm.de/fileadmin/website_uni_ulm/mawi.inst.070/funken/bachelorarbeiten/bachelorthesis_pfh.pdf
    =#

    Ku = size(tabQLeg,1)

    logINVTOL = log(10^(14)) # Logarithm of the inverse of the tolerance of the algorithm

    u, v = real(omg), imag(omg) # Real/Imag part of the frequency
    u2, v2 = u^(2), v^(2) # Squared quantities

    # use an ellipse to decide if we are going up or down in the recursion
    b = min(1.0,4.5/((Ku+1.0)^(1.17))) # Value of ellise parameter
    b2 = b^(2)                          # Squared quantity
    a2 = 1.0+b2                         # Second parameter of the ellipse

    # Testing whether or not we are close to the axis [-1,1]
    if ((u2/a2) + (v2/b2) <= 1.0)       # Within the ellipse
        tabLeg_UP!(omg,val_0,val_1,Ku,tabQLeg) # If we are sufficiently close to [-1,1], we use the forward recurrence
    else                                # Outside of the ellipse
        # Determining the size of the warm-up region
        z2 = abs(u) + abs(v)*im
        ctmp0 = sqrt(z2*z2 - 1.0)
        tmp00 = 2.0*log(abs(z2 + ctmp0))

        K_c = convert(Int64,Ku + ceil(logINVTOL/tmp00)) # Index at which the backward warm-up is started

        tabLeg_BACK!(omg,val_0,K_c,Ku,tabQLeg) # If we are sufficiently far away from [-1,1], we use the backward recurrence
    end
end

function tabPLeg!(omg::Complex{Float64},
                  val_0::Complex{Float64},
                  val_1::Complex{Float64},
                  Ku::Int64,
                  tabPLeg::Array{Complex{Float64},1})
    #=tabPLeg
     Function that pre-computes the Legendre functions, P_k(w),
     for a given complex frequency

     Arguments:
     + omg: COMPLEX frequency. ATTENTION, has to be complex.
     + val_0: Initial value in k=0. ATTENTION, has to be complex.
     + val_1: Initial value in k=1. ATTENTION, has to be complex.
     + tabPLeg: Container where to store the results
     For the Legendre functions, we always use the forward recurrence
    =#
    #####
    # For P_k(w), we always use the UPWARD recurrence.
    tabLeg_UP!(omg,val_0,val_1,Ku,tabPLeg)
end

function tabLeg_UP!(omg::Complex{Float64},
                    val_0::Complex{Float64},
                    val_1::Complex{Float64},
                    Ku::Int64,
                    tabLeg::Array{Complex{Float64},1})
    #=
     Function to compute Legendre functions
     with an UPWARD recurrence
     This function is used for both P_k(w) and Q_k(w)
     as they satisfy the same two-term recurrence
     + It is always used for P_k(w)
     + It is used for Q_k(w) only sufficiently close to [-1,1]

     Arguments:
     + omg: COMPLEX frequency. ATTENTION, has to be complex.
     + val_0: Initial value for k=0. ATTENTION, has to be complex.
     + val_1: Initial value for k=1. ATTENTION, has to be complex.
     + tabLeg: Container where to store the result
    =#
    tabLeg[0+1] = val_0 # Filling in the value of D_0(omg). ATTENTION, to the shift in the array index
    tabLeg[1+1] = val_1 # Filling in the value of D_1(omg). ATTENTION, to the shift in the array index
    #####
    v0, v1 = val_0, val_1 # Initialising the temporary variables used in the recurrence
    #####
    for k=2:(Ku-1) # Loop over the index 2 <= k
        v = ((2.0*k-1.0)*omg*v1 - (k-1.0)*v0)/(k) # Bonnet's recc. relation to get D_k(omg)
        #####
        tabLeg[k+1] = v # Filling in the value of D_k(omg). ATTENTION, to the shift in the array index
        #####
        v0, v1 = v1, v # Shifting the temporary variables
    end
end


function tabLeg_BACK!(omg::Complex{Float64},
                      val_0::Complex{Float64},
                      K_c::Int64,
                      Ku::Int64,
                      tabLeg::Array{Complex{Float64},1})
    #=tabLeg_BACK!
     Function to compute Legendre functions
     with a BACKWARD recurrence

     This function is only used for Q_k(w)
     sufficiently far away from [-1,1].

     Arguments:
     + omg: COMPLEX frequency. ATTENTION, has to be complex.
     + val_0: Initial value in k=0. ATTENTION, has to complex.
     + K_c: Value of k at which the warm-up starts
     + tabLeg: Container where to store the results (modify in place)
    =#

    # Initialisation of D_(K_c+1)[omg]
    v0 = 1.0

    # Initialisation of D_(K_c+2)[omg]
    v1 = 0.0

    # Warm-up phase
    for k=K_c:-1:Ku # ATTENTION, the step is `-1'.
        v = ((2.0*k+3.0)*omg*v0 - (k+2.0)*v1)/(k+1.0) # Computing D_(k)[omg]
        v0, v1 = v, v0 # Shifting the temporary variables
    end

    # Starting to store the values of D_(k)[omg]
    for k=(Ku-1):-1:0 # ATTENTION, the step is `-1'.
        v = ((2.0*k+3.0)*omg*v0 - (k+2.0)*v1)/(k+1.0) # Computing D_(k)[omg]

        tabLeg[k+1] = v # Filling in the value of D_(k)[omg]. ATTENTION, to the shift in the array index

        v0, v1 = v, v0 # Shifting the temporary variables
    end

    pref = val_0/tabLeg[0+1] # Prefactor by which all the values have to be rescaled, so that it is equal to val_0 for k=0. ATTENTION, to the shift of the array.

    for k=0:(Ku-1) # Rescaling all the computed values
        tabLeg[k+1] *= pref # Performing the rescaling. ATTENTION, to the shift in the array index
    end
end

# bring in the individual cases
include("Unstable.jl")
include("Neutral.jl")
include("Damped.jl")


function GetaXi(FHT::structLegendreFHTtype,
                tabGXi::Array{Float64})

    # start with no warnings
    warnflag = zeros(FHT.Ku)

    # start with zero contribution
    res = zeros(FHT.Ku)

    # Loop over the Legendre functions
    for k=1:FHT.Ku

        for i=1:FHT.Ku # Loop over the G-L nodes

            Gval = tabGXi[i] # Current value of G[u_i]

            # check for NaN contribution: skip this contribution in that case
            if isnan(Gval)
                warnflag[k] += 1
                continue
            end

            # check for INF contribution: skip the contribution in that case
            if isinf(Gval)
                warnflag[k] += 1
                continue
            end

            # Current weight
            w = FHT.tabw[i]

            # Current value of P_k
            P = FHT.tabP[k,i]

            res[k] += w*G*P # Update of the sum
        end

        res[k] *= FHT.tabc[k] # Multiplying by the Legendre prefactor.

    end

    return res,warnflag

end
