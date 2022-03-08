#=
 the plasma model
=#

function get_G(u::Float64,qSELF::Float64,xmax::Float64)
    x = u*xmax # Value of the velocity
    return qSELF/(sqrt(pi))*x*exp(-x^(2)) # Returning the value of G(u)
end
