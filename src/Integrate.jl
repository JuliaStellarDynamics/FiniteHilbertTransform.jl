"""
    compute_tabG(uNodes, Gfun)

Compute the G function values on the u nodes.
"""
function compute_tabG(uNodes::Vector{Float64},Gfun::Function)

    K_u = size(uNodes,1)        # Nodes number
    tabG = zeros(Float64,K_u)   # Values array

    Threads.@threads for i=1:K_u # Loop over the nodes

        # Current node position
        u_i = uNodes[i]

        # Compute the value of G[u_i]
        tabG[i] = Gfun(u_i)
    end

    return tabG
end
