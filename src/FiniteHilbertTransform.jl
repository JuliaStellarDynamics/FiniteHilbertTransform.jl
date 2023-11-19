"""
# Finite Hilbert Transform

The Finite Hilbert Transform (FHT) is a discrete mathematical operation that transforms a given function `f(x)` into another function `g(x)` within a finite interval. Unlike the classical Hilbert Transform, which operates on continuous functions over the entire real line, the FHT deals with discrete data points within a specific finite range.

## Applications

- **Signal Processing**: Analyzing and manipulating signals and waveforms.
- **Physics**: Studying wave propagation and interference phenomena.
- **Engineering**: Analyzing electrical circuits, communication systems, and control systems.
"""

"""

For a given discrete function f(x) defined on a finite interval [-1, 1], the FHT transforms it into g(x) using a set of weights and sampling points. Different families of orthogonal functions, such as Legendre polynomials, Chebyshev polynomials, or other orthogonal polynomials, can be used to define the FHT.


"""
module FiniteHilbertTransform



# make an abstract FiniteHilbertTransform type
# @WARNING: Should be defined before any basis definition
abstract type AbstractFHT  end

# bring in the generic integration tools
include("Integrate.jl")

include("IO.jl")
#export dump_tabIminusXi



end # module
