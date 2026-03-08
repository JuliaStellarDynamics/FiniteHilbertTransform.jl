

"""
    get_Xi_DAMPED(n::Int, alpha::Real)

Compute the damping factor for the Chebyshev polynomial of the first kind.

# Arguments
- `n::Int`: The degree of the Chebyshev polynomial.
- `alpha::Real`: The damping parameter.

# Returns
- `Xi::Vector{Complex{Float64}}`: The damping factor for the Chebyshev polynomial.

# Examples
"""
function get_Xi_DAMPED(omg::ComplexF64,struct_tabCheb::ChebyshevFHT)
    sumT, sumU = get_sumT(omg,struct_tabCheb.taba), get_sumU(omg,struct_tabCheb.taba) # Computing the needed sum
    #####
    Xi = -sumT # Starting to compute the expression
    #####
    romg = real(omg) # Real part of the frequency
    #####
    if (romg < -1.0) # On the left of the interval
        Xi -= im*sqrt(1.0-omg^(2))*sumU
    elseif (-1.0 <= romg <= 1.0) # On the right of the interval
        Xi += im*sqrt(1.0-omg^(2))*sumU
    else
        Xi -= im*sqrt(1.0-omg^(2))*sumU
    end
    #####
    return Xi # Output
end


"""
    get_Xi_array(omg::ComplexF64, taba::Vector{Float64})

Compute the Xi array for a given complex frequency `omg` and a vector `taba`.

# Arguments
- `omg::ComplexF64`: The complex frequency.
- `taba::Vector{Float64}`: The input vector.

# Returns
- `Xi`: The computed Xi array.

# Description
This function computes the Xi array by performing a series of calculations based on the given complex frequency `omg` and the input vector `taba`. It first calculates the sumT and sumU values using the `get_sumT` and `get_sumU` functions, respectively. Then, it proceeds to compute the expression for Xi by subtracting the sumT value. The real part of the frequency is stored in the `romg` variable. Depending on the value of `romg`, different calculations are performed to update the value of Xi. Finally, the computed Xi array is returned as the output.

"""
function get_Xi_array(omg::ComplexF64,
                       taba::Vector{Float64})
    sumT, sumU = get_sumT(omg,taba), get_sumU(omg,taba) # Computing the needed sum
    #####
    Xi = -sumT # Starting to compute the expression
    #####
    romg = real(omg) # Real part of the frequency
    #####
    if (romg < -1.0) # On the left of the interval
        Xi -= im*sqrt(1.0-omg^(2))*sumU
    elseif (-1.0 <= romg <= 1.0) # On the right of the interval
        Xi += im*sqrt(1.0-omg^(2))*sumU
    else
        Xi -= im*sqrt(1.0-omg^(2))*sumU
    end
    #####
    return Xi # Output
end
