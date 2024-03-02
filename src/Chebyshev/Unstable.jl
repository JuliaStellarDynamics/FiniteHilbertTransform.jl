##################################################
# Generic function that computes Xi(omg)
# the response matrix in the UNSTABLE case
##################################################
function get_Xi_UNSTABLE(omg::ComplexF64,struct_tabCheb::ChebyshevFHT)
    sumT, sumU = get_sumT(omg,struct_tabCheb.taba), get_sumU(omg,struct_tabCheb.taba) # Computing the needed sum
    #####
    Xi = -sumT + im*sqrt(1.0-omg^(2))*sumU # Computing the unstable expression
    #####
    return Xi
end
