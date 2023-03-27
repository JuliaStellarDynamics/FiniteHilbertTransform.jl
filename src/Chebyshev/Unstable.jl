##################################################
# Generic function that computes Xi(omg)
# the response matrix in the UNSTABLE case
##################################################
function get_Xi_UNSTABLE(omg::ComplexF64,
                         taba::Vector{Float64})
    sumT, sumU = get_sumT(omg,taba), get_sumU(omg,taba) # Computing the needed sum
    #####
    Xi = -sumT + im*sqrt(1.0-omg^(2))*sumU # Computing the unstable expression
    #####
    return Xi # Output
end
