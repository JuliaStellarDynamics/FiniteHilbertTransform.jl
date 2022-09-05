module FiniteHilbertTransform

# bring in the generic integration tools
include("Integrate.jl")

include("IO.jl")
#export dump_tabIminusXi


# make a generic Basis data type
FHTtype = Union{structLegendreFHTtype,structChebyshevFHTtype}



end # module
