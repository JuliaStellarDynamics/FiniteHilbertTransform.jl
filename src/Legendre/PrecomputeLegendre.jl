"""PrecomputeLegendre.jl

Precompute several items for Legendre integration.

"""


"""
Function to initialise the nodes (u) and weights (u)
of the Gauss-Legendre quadrature
"""
function tabuwGLquad(K_u::Int64)

    # pull from FastGaussQuadrature
    # compute the nodes (u) and weights (w) of the G-L quadrature
    tabuGLquad, tabwGLquad = gausslegendre(K_u)

    return tabuGLquad,tabwGLquad
end

"""
Function to initialise the normalisation constant
INVc = 1/c_k = (2k+1)/2
ATTENTION, it corresponds to the INVERSE.
"""
function tabINVcGLquad(K_u::Int64)

    tabINVcGLquad = zeros(Float64,K_u)

    # loop over the Legendre index
    for k=0:(K_u-1)

        # filling in the value. ATTENTION, to the shift in the array index
        tabINVcGLquad[k+1] = (2.0*k+1.0)/(2.0)
    end
    
    return tabINVcGLquad
end

"""
Function to initialise the values of the Legendre polynomials,
as stored in tabPGLquad = P[k,i] = P_k(u_i)
"""
function tabPGLquad(K_u::Int64,tabuGLquad::Vector{Float64})

    tabPGLquad = zeros(Float64,K_u,K_u)
    for i=1:K_u # Loop over the nodes
        u = tabuGLquad[i] # Current value of the node
        #####
        v0 = 1.0 # Initialisation of P_0(u)
        v1 = u   # Initialisation of P_1(u)
        #####
        tabPGLquad[1,i] = v0 # Filling in the value of P_0(u). ATTENTION, to the shift in the array index.
        tabPGLquad[2,i] = v1 # Filling in the value of P_1(u). ATTENTION, to the shift in the array index.
        #####
        for k=2:(K_u-1) # Loop over the index 2 <= k
            v = (2.0*k-1.0)/(k)*u*v1 - (k-1.0)/(k)*v0 # Bonnet's recc. relation to get P_k(u)
            #####
            tabPGLquad[k+1,i] = v # Filling in the value P_k(u). ATTENTION, to the shift in the array index for k
            #####
            v0, v1 = v1, v # Shifting the temporary variables
        end
    end
    return tabPGLquad
end

"""tabGLquad
initialise the nodes and weights of Gauss-Legendre quadrature
"""
function tabGLquad(K_u::Int64)

    # Computing the nodes and weights of the G-L quadrature.
    tabuGLquad,tabwGLquad = tabuwGLquad(K_u)

    # Computing the inverse normalisation of the Legendre polynomials. ATTENTION, corresponds to the inverse of the coefficients.
    INVcGLquad = tabINVcGLquad(K_u)

    # Computing the values of the Legendre polynomials used in the G-L quadrature.
    PGLquad    = tabPGLquad(K_u,tabuGLquad)

    return tabuGLquad,tabwGLquad,INVcGLquad,PGLquad
end
