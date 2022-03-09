

function get_tabuCheby(K_u::Int64)
    tabuCheby = zeros(Float64,K_u)

    for i=0:(K_u-1) # Loop over the nodes. ATTENTION, here we start the loop in 0
        tabuCheby[i+1] = sin((pi)/(2.0*K_u)*(K_u - 2*i - 1)) # Filling in the element. ATTENTION, the array starts in 1
    end

    return tabuCheby
end
