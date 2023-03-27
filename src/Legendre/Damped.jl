


"""
tabLeg!_DAMPED
 Function to fill in all the Legendre arrays
 for a DAMPED mode, i.e. Im[w] < 0.0
"""
function tabLeg!_DAMPED(omg::ComplexF64,
                        struct_tabLeg::LegendreFHT)

    # the container for D_k(w)
    tabDLeg = struct_tabLeg.tabDLeg

    # the container for Q_k(w)
    tabQLeg = struct_tabLeg.tabQLeg

    # now, compute Q_k(w)

    # fill the initial values of Q_k(w)
    val_0_Q = log(1.0-omg) - log(-1.0-omg)
    val_1_Q = 2.0 + omg*val_0_Q

    # compute the rest of the k values
    tabQLeg!(omg,val_0_Q,val_1_Q,tabQLeg)

    # compute the value of H[-1,1,Re[w]].
    val_H = Heaviside(real(omg))

    # if we are outside of the interval [-1,1], we do not need to compute P_k(w) (it =0)
    if (val_H == 0.0)

        # Computing D_k(w)
        for k=0:(struct_tabLeg.Ku-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays.
        end

    # if we are within the interval [-1,1], we need to compute P_k(w)
    else

        # the container for P_k(w)
        tabPLeg = struct_tabLeg.tabPLeg

        # Initial values of P_k(w)
        val_0_P = 1.0  + 0.0*im # Initial value for k = 0. ATTENTION, we make sure that the quantity is seen as complex.
        val_1_P = omg           # Initial value for k = 1.
        #####
        tabPLeg!(omg,val_0_P,val_1_P,struct_tabLeg.Ku,tabPLeg) # Computing the tabPLeg
        ##########
        # Computing the D_k(w)
        #####
        for k=0:(struct_tabLeg.Ku-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] + 2.0*im*pi*val_H*tabPLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays
        end
        #####
    end
end
