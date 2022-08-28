


# Bring in the prefactors
#include("Legendre/PrecomputeLegendre.jl")

# Bring in the integration tools
include("Rational/Rational.jl")



function tabIminusXi!()
    Threads.@threads for iomega=1:nomega # Loop over all the considered COMPLEX frequencies
        #####
        omg = tabomega[iomega] # Current value of the complex frequency
        #####
        taba = struct_Rational_serial.taba # Coefficients of the continued fraction to use in the rational function
        #####
        val = evaluate_Rational(omg,tabomega_Rational,taba) # Computing I-Xi(omg)
        #####
        tabIminusXi[iomega] = val # Filling in tabIminusXi
    end
end
