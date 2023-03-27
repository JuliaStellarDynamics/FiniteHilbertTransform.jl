module FiniteHilbertTransform

# make an abstract FiniteHilbertTransform type
# @WARNING: Should be defined before any basis definition
abstract type AbstractFHT  end

# bring in the generic integration tools
include("Integrate.jl")

include("IO.jl")
#export dump_tabIminusXi



end # module
