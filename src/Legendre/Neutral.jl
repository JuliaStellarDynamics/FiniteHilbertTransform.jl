##################################################
# Function to fill in all the Legendre arrays
# for an NEUTRAL mode, i.e. Im[w] = 0.0
##################################################
function tabLeg!_NEUTRAL(omg::Complex{Float64},
                         struct_tabLeg::structLegendreFHTtype)
    ##########
    tabDLeg = struct_tabLeg.tabDLeg # Name of the container for D_k(w)
    ##########
    romg = real(omg) # Keeping only the real part of omg
    ##########
    # Computing Q_k(w)
    #####
    tabQLeg = struct_tabLeg.tabQLeg # Name of the container for Q_k(w)
    ######
    # Initial values of Q_k(w)
    # The integral is computed in the Cauchy sense, hence the absolute values.
    val_0_Q = log(abs(1.0-romg)) - log(abs(-1.0-romg)) + 0.0*im # Initial value for k = 0. ATTENTION, we make sure that the quantity is seen as complex. ATTENTION, not to forget the `abs'.
    val_1_Q = 2.0 + romg*val_0_Q                       + 0.0*im # Initial value for k = 1. ATTENTION, we make sure that the quantity is seen as complex.
    #####
    tabQLeg!(omg,val_0_Q,val_1_Q,tabQLeg) # Computing the tabQLeg
    ##########
    # Computing P_k(w) if needed
    #####
    val_H = Heaviside(romg) # Computing the value of H[-1,1,Re[w]]. ATTENTION, the argument must be real.
    #####
    if (val_H == 0.0) # We are outside of the interval [-1,1]: no need to compute P_k(w)
        ##########
        # Computing D_k(w)
        ##########
        for k=0:(struct_tabLeg.K_u-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays.
        end
        #####
    else # We are within the interval [-1,1]: we need to compute P_k(w)
        #####
        tabPLeg = struct_tabLeg.tabPLeg # Name of the container for P_k(w)
        #####
        # Initial values of P_k(w)
        val_0_P = 1.0  + 0.0*im # Initial value for k = 0. ATTENTION, we make sure that the quantity is seen as complex.
        val_1_P = romg + 0.0*im # Initial value for k = 1. ATTENTION, we make sure that the quantity is seen as complex.
        #####
        tabPLeg!(omg,val_0_P,val_1_P,struct_tabLeg.K_u,tabPLeg) # Computing the tabPLeg
        ##########
        # Computing D_k(w)
        #####
        for k=0:(struct_tabLeg.K_u-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] + im*pi*val_H*tabPLeg[k+1] # Filling in tabDLeg. ATTENTION, to the shift of the arrays
        end
        #####
    end
end
