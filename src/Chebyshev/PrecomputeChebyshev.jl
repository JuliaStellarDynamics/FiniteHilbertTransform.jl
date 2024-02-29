



function get_tabuCheby(Ku::Int64)
    tabuCheby = zeros(Float64,Ku)

    for i=0:(Ku-1) # Loop over the nodes. ATTENTION, here we start the loop in 0
        tabuCheby[i+1] = sin((pi)/(2.0*Ku)*(Ku - 2*i - 1)) # Filling in the element. ATTENTION, the array starts in 1
    end

    return tabuCheby
end


"""

sampling points for Chebyshev integration
"""
function GettabuCquad(Ku::Int64)
    tabuCheby = zeros(Float64,Ku)

    for i=0:(Ku-1) # Loop over the nodes. ATTENTION, here we start the loop in 0
        tabuCheby[i+1] = sin((pi)/(2.0*Ku)*(Ku - 2*i - 1)) # Filling in the element. ATTENTION, the array starts in 1
    end

    return tabuCheby
end




"""tabCquad
initialise the nodes and weights of Chebyshev quadrature
"""
function tabCquad(Ku::Int64)

    # Computing the nodes and weights of the Chebyshev quadrature.
    tabuCquad = GettabuCquad(Ku)

    tabwCquad = ones(Ku)

    # only extra normalisation for Chebyshev comes from taking the FFT
    INVcCquad = ones(Ku) ./ Ku

    # Computing the values of the Chebyshev polynomials used in the Chebyshev quadrature.
    PCquad    = getTu(Ku,tabuCquad)

    return tabuCquad,tabwCquad,INVcCquad,PCquad
end
