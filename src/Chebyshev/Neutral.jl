##################################################
# Response matrix in the NEUTRAL case
# ATTENTION, we do not enforce for omg
# to have a vanishing imaginary part
##################################################
function get_Xi_NEUTRAL(omg::ComplexF64,
                         taba::Vector{Float64})
    sumT, sumU = get_sumT(omg,taba), get_sumU(omg,taba) # Computing the needed sum
    #####
    Xi = -sumT # Starting to compute the expression
    #####
    romg = real(omg) # Real part of the frequency
    #####
    if (romg < -1.0) # On the left of the interval
        Xi -=    sqrt(omg^(2)-1.0)*sumU
    elseif (-1.0 <= romg <= 1.0) # On the right of the interval
        Xi += im*sqrt(1.0-omg^(2))*sumU
    else
        Xi +=    sqrt(omg^(2)-1.0)*sumU
    end
    #####
    return Xi # Output
end
