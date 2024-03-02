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
    LegendreFHT

Type representing Legendre Finite Hilbert Transform (FHT) parameters.

This struct stores the necessary parameters for performing a Finite Hilbert Transform using Legendre functions. It implements the AbstractFHT interface and provides the required data structures and computations for Legendre FHT.

# Fields:
- `name::String`: FHT name (default: "Legendre").
- `Ku::Int64`: Number of sample points.
- `tabu::Array{Float64,1}`: Array containing u values (sampling points).
- `tabw::Array{Float64,1}`: Array containing w values (weights at sampling points).
- `tabP::Matrix{Float64}`: Matrix containing P_k(u) values (Ku x Ku).
- `tabc::Vector{Float64}`: Vector containing prefactor values at each sampling point.
- `tabPLeg::Array{ComplexF64,1}`: Static container for tabPLeg (Legendre functions of the first kind).
- `tabQLeg::Array{ComplexF64,1}`: Static container for tabQLeg (Hilbert-transformed Legendre functions).
- `tabDLeg::Array{ComplexF64,1}`: Static container for tabDLeg (Derivatives of Legendre functions).

# Example:
```julia
Ku = 100
tabu = collect(-1.0:2/(Ku-1):1.0)
tabw = compute_weights(tabu) # Compute weights for Legendre functions
tabP = compute_legendre_matrix(tabu) # Compute Legendre functions matrix
tabc = compute_prefactors(tabu) # Compute Legendre prefactors

FHT = LegendreFHT("Legendre", Ku, tabu, tabw, tabP, tabc, zeros(ComplexF64, Ku), zeros(ComplexF64, Ku), zeros(ComplexF64, Ku))
```
"""
struct LegendreFHT <: AbstractFHT

    name::String         # FHT name (default Legendre)
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
    LegendreFHT(Ku[, name])

Create a LegendreFHT structure.

# Arguments
- `Ku::Int64`: Number of Legendre modes.
- `name::String`: Name of the LegendreFHT structure (default: "Legendre").

# Returns
- `legendre_fht::LegendreFHT`: LegendreFHT structure.

# Description
`LegendreFHT` constructs a LegendreFHT structure with the specified number of Legendre modes (`Ku`). Optionally, you can provide a `name` for the structure.

"""
function LegendreFHT(Ku::Int64;name::String="Legendre")

    tabu,tabw,tabc,tabP = tabGLquad(Ku)

    return LegendreFHT(name,Ku,tabu,tabw,tabP,tabc,zeros(ComplexF64,Ku),zeros(ComplexF64,Ku),zeros(ComplexF64,Ku))

end



"""
    GettabD!(omg::Complex{Float64}, struct_tabLeg::LegendreFHT; verbose::Int64=0)

Fill `struct_tabLeg` at a given complex frequency `omg` for the integration being considered.

# Arguments
- `omg::Complex{Float64}`: Complex frequency value.
- `struct_tabLeg::LegendreFHT`: LegendreFHT structure to fill.
- `verbose::Int64`: Verbosity level (default: 0).

# Description
`GettabD!` automatically selects the integration style based on the imaginary part of `omg`. If the imaginary part is negative, it uses damped integration. If it's exactly zero, it uses neutral mode calculation. Otherwise, it uses unstable integration.

# Example
```julia
tabLeg = LegendreFHT(10)
omg = 1.0 + 0.5im
GettabD!(omg, tabLeg, verbose=2)
```

"""
function GettabD!(omg::ComplexF64,
                  struct_tabLeg::LegendreFHT;
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


"""
    Heaviside(x::Float64)

Calculate the Heaviside function on the interval [-1, 1].

The Heaviside function, denoted as H(x), is defined as follows:
- H(x) = 0 for x < -1
- H(x) = 0.5 for x = -1
- H(x) = 1 for -1 < x < 1
- H(x) = 0.5 for x = 1
- H(x) = 0 for x > 1

ATTENTION: The equality tests on Float64 might not be very robust.

# Arguments:
- `x::Float64`: Real number for which Heaviside function is calculated.

# Returns:
- `Float64`: Value of the Heaviside function at x.

# Example:
```julia
x = 0.5
result = Heaviside(x)  # Returns 1.0
```
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

"""
    tabQLeg!(omg::ComplexF64, val_0::ComplexF64, val_1::ComplexF64, tabQLeg::Array{ComplexF64,1})

Precompute Hilbert-transformed Legendre functions for a given complex frequency.


# Arguments:
- `omg::ComplexF64`: Complex frequency. (ATTENTION: Must be complex.)
- `val_0::ComplexF64`: Initial value for ` k = 0 `. (ATTENTION: Must be complex.)
- `val_1::ComplexF64`: Initial value for ` k = 1 `. (ATTENTION: Must be complex.)
- `tabQLeg::Array{ComplexF64,1}`: Container to store the precomputed Hilbert-transformed Legendre functions.

# Details:
This function uses different recurrence relations based on the location of the complex frequency ` omg `. If ` omg ` is sufficiently close to the real line `[-1,1]`, it employs an upward recurrence. Otherwise, if ` omg ` is far away from the real line `[-1,1]`, it uses a backward recurrence. The transition between these regimes is determined dynamically.

The transition from the two regimes follows from the thesis Stable Implementation of Three-Term Recurrence Relations, Pascal Frederik Heiter, June, 2010
https://www.uni-ulm.de/fileadmin/website_uni_ulm/mawi.inst.070/funken/bachelorarbeiten/bachelorthesis_pfh.pdf

# Example:
```julia
omg = 1.0 + 2.0im
val_0 = 1.0 + 1.0im
val_1 = 2.0 + 2.0im
Ku = 10
tabQLeg = zeros(ComplexF64, Ku)
tabQLeg!(omg, val_0, val_1, tabQLeg)
```

"""
function tabQLeg!(omg::ComplexF64,
                  val_0::ComplexF64,
                  val_1::ComplexF64,
                  tabQLeg::Array{ComplexF64,1})

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

"""
    tabPLeg!(omg::Complex{Float64}, val_0::Complex{Float64}, val_1::Complex{Float64}, Ku::Int64, tabPLeg::Vector{Complex{Float64}})

Pre-computes the Legendre functions, P_k(w), for a given complex frequency.

# Arguments
- `omg::Complex{Float64}`: Complex frequency.
- `val_0::Complex{Float64}`: Initial value in k=0.
- `val_1::Complex{Float64}`: Initial value in k=1.
- `Ku::Int64`: Number of Legendre modes.
- `tabPLeg::Vector{Complex{Float64}}`: Container to store the results.

# Description
`tabPLeg!` pre-computes the Legendre functions, P_k(w), for a given complex frequency `omg` using the upward recurrence method. The results are stored in `tabPLeg`.

# Example
```julia
omg = 1.0 + 0.5im
val_0 = 0.0 + 0.0im
val_1 = 1.0 + 0.0im
Ku = 10
tabPLeg = zeros(Complex{Float64}, Ku)
tabPLeg!(omg, val_0, val_1, Ku, tabPLeg)
```

"""
function tabPLeg!(omg::ComplexF64,
                  val_0::ComplexF64,
                  val_1::ComplexF64,
                  Ku::Int64,
                  tabPLeg::Array{ComplexF64,1})

    tabLeg_UP!(omg,val_0,val_1,Ku,tabPLeg)
end


"""
    tabLeg_UP!(omg::ComplexF64, val_0::ComplexF64, val_1::ComplexF64, Ku::Int64, tabLeg::Array{ComplexF64,1})

Compute Legendre functions with an UPWARD recurrence.

This function calculates Legendre functions with an upward recurrence relation. Both P_k(w) and Q_k(w) satisfy the same two-term recurrence, and this function is used for both P_k(w) and Q_k(w) computations.

# Arguments:
- `omg::ComplexF64`: Complex frequency. (ATTENTION: Must be complex.)
- `val_0::ComplexF64`: Initial value for k=0. (ATTENTION: Must be complex.)
- `val_1::ComplexF64`: Initial value for k=1. (ATTENTION: Must be complex.)
- `Ku::Int64`: Upper limit of the Legendre functions to be computed.
- `tabLeg::Array{ComplexF64,1}`: Container to store the resulting Legendre functions.

# Details:
This function computes Legendre functions using Bonnet's recurrence relation: D_k(omg) = ((2.0*k-1.0)*omg*v_1 - (k-1.0)*v_0) / k, where v_0 and v_1 are initial values.

# Example:
```julia
omg = 1.0 + 2.0im
val_0 = 1.0 + 1.0im
val_1 = 2.0 + 2.0im
Ku = 10
tabLeg = zeros(ComplexF64, Ku)
tabLeg_UP!(omg, val_0, val_1, Ku, tabLeg)
```
"""
function tabLeg_UP!(omg::ComplexF64,
                    val_0::ComplexF64,
                    val_1::ComplexF64,
                    Ku::Int64,
                    tabLeg::Array{ComplexF64,1})

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


"""
    tabLeg_BACK!(omg::ComplexF64, val_0::ComplexF64, K_c::Int64, Ku::Int64, tabLeg::Array{ComplexF64,1})

Compute Legendre functions with a BACKWARD recurrence.

This function calculates Legendre functions with a backward recurrence relation. It is specifically used for computing Q_k(w) when the input frequency omg is sufficiently far away from the range [-1,1].

# Arguments:
- `omg::ComplexF64`: Complex frequency. (ATTENTION: Must be complex.)
- `val_0::ComplexF64`: Initial value for k=0. (ATTENTION: Must be complex.)
- `K_c::Int64`: Value of k at which the warm-up starts.
- `Ku::Int64`: Upper limit of the Legendre functions to be computed.
- `tabLeg::Array{ComplexF64,1}`: Container to store the resulting Legendre functions. (Modified in place.)

# Details:
This function computes Legendre functions using a backward recurrence relation. It performs a warm-up phase for k in the range K_c down to Ku, and then stores the computed Legendre values from Ku down to 0. The computed values are rescaled to match the initial value val_0 for k=0.

# Example:
```julia
omg = 1.0 + 2.0im
val_0 = 1.0 + 1.0im
K_c = 5
Ku = 10
tabLeg = zeros(ComplexF64, Ku)
tabLeg_BACK!(omg, val_0, K_c, Ku, tabLeg)
```
"""
function tabLeg_BACK!(omg::ComplexF64,
                      val_0::ComplexF64,
                      K_c::Int64,
                      Ku::Int64,
                      tabLeg::Array{ComplexF64,1})


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


"""
    GetaXi!(FHT::LegendreFHT, tabGXi::AbstractVector{Float64}, res::Vector{Float64}, warnflag::Vector{Float64})

Compute the Finite Hilbert Transform for Legendre functions.

# Arguments
- `FHT::LegendreFHT`: An object representing the Legendre Finite Hilbert Transform.
- `tabGXi::AbstractVector{Float64}`: Vector containing precomputed values of G[u_i].
- `res::Vector{Float64}`: Output vector to store the results of the Finite Hilbert Transform.
- `warnflag::Vector{Float64}`: Vector to store warning flags for each Legendre function.

# Details
This function computes the Finite Hilbert Transform for Legendre functions based on the provided precomputed values of G[u_i].

# Output
- `res::Vector{Float64}`: Vector containing the results of the Finite Hilbert Transform for Legendre functions.
- `warnflag::Vector{Float64}`: Vector containing warning flags. Each element represents the number of NaN or INF contributions for the corresponding Legendre function.

# Example
```julia
FHT = LegendreFHT(parameters)
tabGXi = compute_tabGXi(...)  # Precompute tabGXi values
res = zeros(Float64, FHT.Ku)
warnflag = zeros(Float64, FHT.Ku)
GetaXi!(FHT, tabGXi, res, warnflag)
```
"""
function GetaXi!(FHT::LegendreFHT,
                 tabG::AbstractVector{Float64},
                 res::Vector{Float64},warnflag::Int64)

    # check vector length for agreement before proceeding
    @assert length(tabG) == FHT.Ku "FiniteHilbertTransform.Legendre.GetaXi!: tabG length is not the same as Ku."

    # check vector for nan or inf: if there are, loop through and zero out.
    # see speed discussion here: we might be able to do even better
    # https://discourse.julialang.org/t/fastest-way-to-check-for-inf-or-nan-in-an-array/76954/27
    if 1==0#isfinite(sum(tabG))

        # Loop over the Legendre functions
        for k=1:FHT.Ku

            # do the inner sum over u as a vector

            # weight * G(u) * P_k(u)
            res[k] = sum(FHT.tabw .* tabG .* FHT.tabP[k]) # Update of the sum

            res[k] *= FHT.tabc[k] # Multiplying by the Legendre prefactor.

        end

    else

        # Loop over the Legendre functions
        for k=1:FHT.Ku

            res[k] = 0.0

            for i=1:FHT.Ku # Loop over the G-L nodes

                Gval = tabG[i] # Current value of G[u_i]

                # check for NaN contribution: skip this contribution in that case
                if isnan(Gval) || isinf(Gval)
                    warnflag += 1
                    continue
                end

                # Current weight
                w = FHT.tabw[i]

                # Current value of P_k
                P = FHT.tabP[k,i]

                res[k] += w*Gval*P # Update of the sum
            end

            res[k] *= FHT.tabc[k] # Multiplying by the Legendre prefactor.

        end

    end

    return res,warnflag

end



"""
    GetaXi(FHT::LegendreFHT, tabGXi::Array{Float64})

Calculate the Finite Hilbert Transform for Legendre functions.

This function computes the Finite Hilbert Transform for Legendre functions based on the provided Legendre Finite Hilbert Transform parameters and precomputed values of G[u_i]. The results are stored in a vector, and warning flags indicating NaN or INF contributions are also provided.

# Arguments:
- `FHT::LegendreFHT`: An object representing the Legendre Finite Hilbert Transform.
- `tabGXi::Array{Float64}`: Array containing precomputed values of G[u_i].

# Returns:
- `res::Vector{Float64}`: Vector containing the results of the Finite Hilbert Transform for Legendre functions.
- `warnflag::Vector{Float64}`: Vector containing warning flags. Each element represents the number of NaN or INF contributions for the corresponding Legendre function.

# Example:
```julia
FHT = LegendreFHT(parameters)
tabGXi = compute_tabGXi(...)  # Precompute tabGXi values
res, warnflag = GetaXi(FHT, tabGXi)
```
"""
function GetaXi(FHT::LegendreFHT,
                tabGXi::Array{Float64})

    # start with no warnings
    warnflag = 0

    # start with zero contribution
    res = zeros(FHT.Ku)

    GetaXi!(FHT,tabGXi,res,warnflag)

    return res,warnflag
end


"""
Compute the values of I-Xi(omg) for a given complex frequency.

# Arguments
- `omg::Complex{Float64}`: Complex frequency.
- `taba::Vector{Float64}`: Vector of coefficients a_k(u).
- `xmax::Float64`: Maximum value of x.
- `FHT::FiniteHilbertTransform.AbstractFHT`: AbstractFHT structure.

# Returns
- `IminusXi::Complex{Float64}`: Value of I-Xi(omg).

# Example
```julia
GetLegendreIminusXiPlasma(1.0 + 1.0im, taba, 10.0, FHT)
```
"""
function GetIminusXi(ϖ::Complex{Float64},taba::Vector{Float64},FHT::LegendreFHT)


    # Rescale the COMPLEX frequency
    K_u = size(taba,1)

    # compute the Hilbert-transformed Legendre functions
    FiniteHilbertTransform.GettabD!(ϖ,FHT)

    xi = 0.0 + 0.0*im # Initialise xi

    # loop over the Legendre functions
    for k=1:(K_u)

        # add the contribution
        xi += taba[k]*FHT.tabDLeg[k]
    end

    # compute 1.0 - xi
    IminusXi = 1.0 - xi

    return IminusXi # Output
end

