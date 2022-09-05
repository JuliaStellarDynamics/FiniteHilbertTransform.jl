# perform discrete sine transform in Gauss-Chebyshev quadrature
using FFTW

# access to the nodes and weights of the Gauss-Legendre quadrature
using FastGaussQuadrature

# Chebyshev:
# Bring in the prefactors
include("Chebyshev/PrecomputeChebyshev.jl")

# Bring in the integration tools
include("Chebyshev/Chebyshev.jl")

# Gauss-Legendre:
# Bring in the prefactors
include("Legendre/PrecomputeLegendre.jl")

# Bring in the integration tools
include("Legendre/Legendre.jl")
