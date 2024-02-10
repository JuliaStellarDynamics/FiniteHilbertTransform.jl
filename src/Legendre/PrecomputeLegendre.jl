"""PrecomputeLegendre.jl

Precompute several items for Legendre integration.

"""



"""
    tabuwGLquad(K_u::Int64)

Initialize the nodes (u) and weights (w) of the Gauss-Legendre quadrature.

# Arguments
- `K_u::Int64`: Number of Legendre modes.

# Returns
- `tabuGLquad::Vector{Float64}`: Vector of nodes (u) of the Gauss-Legendre quadrature.
- `tabwGLquad::Vector{Float64}`: Vector of weights (w) of the Gauss-Legendre quadrature.

# Description
`tabuwGLquad` initializes the nodes (u) and weights (w) of the Gauss-Legendre quadrature using the FastGaussQuadrature package. It computes the nodes and weights for the Gauss-Legendre quadrature and returns them as vectors.

# Example
```julia
K_u = 5
tabuGLquad, tabwGLquad = tabuwGLquad(K_u)
```

# Notes
- This function computes the nodes and weights for the Gauss-Legendre quadrature using external libraries.
"""
function tabuwGLquad(K_u::Int64)
    # Pull from FastGaussQuadrature
    # Compute the nodes (u) and weights (w) of the G-L quadrature
    tabuGLquad, tabwGLquad = gausslegendre(K_u)
    return tabuGLquad, tabwGLquad
end


"""
    tabINVcGLquad(K_u::Int64)

Initialize the normalization constant INVc = 1/c_k = (2k+1)/2.

# Arguments
- `K_u::Int64`: Number of Legendre modes.

# Returns
- `tabINVcGLquad::Vector{Float64}`: Vector of normalization constants.

# Description
`tabINVcGLquad` initializes the normalization constant INVc for Legendre modes. It computes the normalization constant for each Legendre mode and stores the results in a vector.

# Example
```julia
K_u = 5
tabINVcGLquad = tabINVcGLquad(K_u)
```

# Notes
- This function computes the normalization constant for Legendre modes indexed from 0 to K_u-1.
- ATTENTION: It corresponds to the INVERSE normalization constant.

"""
function tabINVcGLquad(K_u::Int64)
    tabINVcGLquad = zeros(Float64, K_u)
    # Loop over the Legendre index
    for k in 0:(K_u-1)
        # Filling in the value. ATTENTION to the shift in the array index
        tabINVcGLquad[k+1] = (2.0 * k + 1.0) / 2.0
    end
    return tabINVcGLquad
end


"""
    tabPGLquad(K_u::Int64, tabuGLquad::Vector{Float64})

Initialize the values of the Legendre polynomials.

# Arguments
- `K_u::Int64`: Number of Legendre modes.
- `tabuGLquad::Vector{Float64}`: Vector of nodes (u) of the Gauss-Legendre quadrature.

# Returns
- `tabPGLquad::Matrix{Float64}`: Matrix of Legendre polynomials.

# Description
`tabPGLquad` initializes the values of the Legendre polynomials and stores them in a matrix. Each column of the matrix corresponds to the Legendre polynomials evaluated at a specific node (u) from the Gauss-Legendre quadrature.

# Example
```julia
K_u = 5
tabuGLquad, _ = tabuwGLquad(K_u)
tabPGLquad = tabPGLquad(K_u, tabuGLquad)
```

# Notes
- This function computes the Legendre polynomials using Bonnet's recurrence relation.
- ATTENTION: It corresponds to the INVERSE normalization constant.
"""
function tabPGLquad(K_u::Int64, tabuGLquad::Vector{Float64})
    tabPGLquad = zeros(Float64, K_u, K_u)
    for i in 1:K_u # Loop over the nodes
        u = tabuGLquad[i] # Current value of the node
        v0 = 1.0           # Initialisation of P_0(u)
        v1 = u             # Initialisation of P_1(u)
        tabPGLquad[1, i] = v0 # Filling in the value of P_0(u)
        tabPGLquad[2, i] = v1 # Filling in the value of P_1(u)
        for k in 2:(K_u-1) # Loop over the index 2 <= k
            v = (2.0 * k - 1.0) / k * u * v1 - (k - 1.0) / k * v0 # Bonnet's recc. relation to get P_k(u)
            tabPGLquad[k + 1, i] = v # Filling in the value P_k(u)
            v0, v1 = v1, v             # Shifting the temporary variables
        end
    end
    return tabPGLquad
end


"""
    tabGLquad(K_u::Int64)

Initialize the nodes and weights of Gauss-Legendre quadrature.

# Arguments
- `K_u::Int64`: Number of Legendre modes.

# Returns
- `tabuGLquad::Vector{Float64}`: Vector of nodes (u) of the Gauss-Legendre quadrature.
- `tabwGLquad::Vector{Float64}`: Vector of weights (w) of the Gauss-Legendre quadrature.
- `INVcGLquad::Vector{Float64}`: Vector of normalization constants for Legendre polynomials.
- `PGLquad::Matrix{Float64}`: Matrix of Legendre polynomials evaluated at quadrature points.

# Description
`tabGLquad` initializes the nodes (u) and weights (w) of the Gauss-Legendre quadrature, along with the normalization constants and Legendre polynomials evaluated at quadrature points.

# Example
```julia
K_u = 5
tabuGLquad, tabwGLquad, INVcGLquad, PGLquad = tabGLquad(K_u)
```

"""
function tabGLquad(K_u::Int64)
    # Computing the nodes and weights of the G-L quadrature.
    tabuGLquad, tabwGLquad = tabuwGLquad(K_u)
    # Computing the inverse normalization of the Legendre polynomials
    INVcGLquad = tabINVcGLquad(K_u)
    # Computing the values of the Legendre polynomials used in the G-L quadrature.
    PGLquad = tabPGLquad(K_u, tabuGLquad)
    return tabuGLquad, tabwGLquad, INVcGLquad, PGLquad
end
