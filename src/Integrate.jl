


function compute_tabG(uNodes::Vector{Float64},
                      qSELF::Float64,
                      xmax::Float64)
    #=
     Pre-computes the needed values of G(u) for the specified plasma model
    =#

    K_u = size(uNodes,1)
    tabG = zeros(Float64,K_u)

    # Loop over the nodes
    Threads.@threads for i=1:K_u

        # Current node position
        u_i = uNodes[i]

        # Compute the value of G[u_i]
        tabG[i] = get_G(u_i,qSELF,xmax)

    end

    return tabG
end
