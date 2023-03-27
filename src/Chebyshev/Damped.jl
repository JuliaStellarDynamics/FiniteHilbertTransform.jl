##################################################
# Response matrix in the DAMPED case
##################################################
function get_Xi_DAMPED(omg::ComplexF64,
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
