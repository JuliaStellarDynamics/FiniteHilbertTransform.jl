



function get_tabuCheby(K_u::Int64)
    tabuCheby = zeros(Float64,K_u)

    for i=0:(K_u-1) # Loop over the nodes. ATTENTION, here we start the loop in 0
        tabuCheby[i+1] = sin((pi)/(2.0*K_u)*(K_u - 2*i - 1)) # Filling in the element. ATTENTION, the array starts in 1
    end

    return tabuCheby
end


"""

sampling points for Chebyshev integration
"""
function tabuCquad(K_u::Int64)
    tabuCheby = zeros(Float64,K_u)

    for i=0:(K_u-1) # Loop over the nodes. ATTENTION, here we start the loop in 0
        tabuCheby[i+1] = sin((pi)/(2.0*K_u)*(K_u - 2*i - 1)) # Filling in the element. ATTENTION, the array starts in 1
    end

    return tabuCheby
end




"""tabCquad
initialise the nodes and weights of Chebyshev quadrature
"""
function tabCquad(Ku::Int64)

    # Computing the nodes and weights of the G-L quadrature.
    tabuCquad,tabwCquad = tabuwCquad(Ku)

    # only extra normalisation for Chebyshev comes from taking the FFT
    INVcCquad = ones(Ku) ./ Ku

    # Computing the values of the Legendre polynomials used in the G-L quadrature.
    PCquad    = getTu(Ku,tabuCquad)

    return tabuCquad,tabwCquad,INVcCquad,PCquad
end
