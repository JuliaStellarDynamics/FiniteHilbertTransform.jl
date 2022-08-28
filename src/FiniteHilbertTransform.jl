module FiniteHilbertTransform

# bring in the generic integration tools
include("Integrate.jl")

include("IntegrateLegendre.jl")
#export setup_legendre_integration,compute_tabIminusXi,get_Legendre_IminusXi,test_ninepointsL

include("IntegrateChebyshev.jl")
#export setup_chebyshev_integration,compute_tabIminusXi,get_Chebyshev_IminusXi,test_ninepointsC

include("IO.jl")
#export dump_tabIminusXi

end # module
