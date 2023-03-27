##################################################
# Function to fill in all the Legendre arrays
# for an UNSTABLE mode, i.e. Im[w] > 0.0
##################################################
function tabLeg!_UNSTABLE(omg::ComplexF64,
                          struct_tabLeg::LegendreFHT)
    ##########
    tabDLeg = struct_tabLeg.tabDLeg # Name of the container for D_k(w)
    ##########
    # Computing Q_k(w)
    #####
    tabQLeg = struct_tabLeg.tabQLeg # Name of the container for Q_k(w)
    ######
    # Initial values of Q_k(w)
    val_0_Q = log(1.0-omg) - log(-1.0-omg) # Initial value for k = 0
    val_1_Q = 2.0 + omg*val_0_Q           # Initial value for k = 1
    #####
    tabQLeg!(omg,val_0_Q,val_1_Q,tabQLeg) # Computing the tabQLeg
    ##########
    # Computing D_k(w)
    #####
    for k=0:(struct_tabLeg.Ku-1) # Loop over the Legendre indices
        tabDLeg[k+1] = tabQLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays
    end
end
