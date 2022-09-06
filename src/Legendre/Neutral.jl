

"""
Function to fill in all the Legendre arrays for an NEUTRAL mode, i.e. Im[w] = 0.0
"""
function tabLeg!_NEUTRAL(omg::Complex{Float64},
                         struct_tabLeg::structLegendreFHTtype)

    # The integral is computed in the Cauchy sense, hence the absolute values.

    # simplify naming convention for the container for D_k(w)
    tabDLeg = struct_tabLeg.tabDLeg
    tabQLeg = struct_tabLeg.tabQLeg
    tabPLeg = struct_tabLeg.tabPLeg

    # keep only the real part of omg
    romg = real(omg)

    # initial values of Q_k(w):

    if romg == 1.0
        # unique case, requires special seed values
        val_0_Q = 0.0 + 0.0im
        val_1_Q = 2.0 + 0.0im

    else
        # initial value for k = 0. ATTENTION, not to forget the `abs'.
        val_0_Q = log(abs(1.0-romg)) - log(abs(-1.0-romg)) + 0.0*im

        # initial value for k = 1.
        val_1_Q = 2.0 + romg*val_0_Q                       + 0.0*im
    end

    # now the tabQLeg
    tabQLeg!(omg,val_0_Q,val_1_Q,tabQLeg)

    # compute P_k(w) if needed:
    # Computing the value of H[-1,1,Re[w]]. ATTENTION, the argument must be real.
    val_H = Heaviside(romg)

    # We are outside of the interval [-1,1]: no need to compute P_k(w)
    if (val_H == 0.0)

        # compute D_k(w):

        # Loop over the Legendre indices
        for k=1:(struct_tabLeg.Ku)

            # Filling in tabDLeg
            tabDLeg[k] = tabQLeg[k]
        end

    # We are within the interval [-1,1]: we need to compute P_k(w)
    else

        # Initial values of P_k(w):
        # Initial value for k = 0.
        val_0_P = 1.0  + 0.0*im

        # Initial value for k = 1.
        val_1_P = romg + 0.0*im

        tabPLeg!(omg,val_0_P,val_1_P,struct_tabLeg.Ku,tabPLeg) # Computing the tabPLeg

        # Computing D_k(w):

        # Loop over the Legendre indices
        for k=1:(struct_tabLeg.Ku)

            # Filling in tabDLeg
            tabDLeg[k] = tabQLeg[k] + im*pi*val_H*tabPLeg[k]
        end

    end
end
