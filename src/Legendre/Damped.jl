"""
    tabLeg!_DAMPED(omg::Complex{Float64}, struct_tabLeg::LegendreFHT)

Fill in all the Legendre arrays for a DAMPED mode, i.e., Im[w] < 0.0.

# Arguments
- `omg::Complex{Float64}`: Complex frequency.
- `struct_tabLeg::LegendreFHT`: LegendreFHT structure.

# Description
`tabLeg!_DAMPED` fills in all the Legendre arrays for a DAMPED mode, where Im[w] < 0.0. It computes the Legendre coefficients D_k(w) and, if needed, Legendre polynomials P_k(w) for the given complex frequency.

"""
function tabLeg!_DAMPED(omg::Complex{Float64}, struct_tabLeg::LegendreFHT)
    # The container for D_k(w)
    tabDLeg = struct_tabLeg.tabDLeg
    # The container for Q_k(w)
    tabQLeg = struct_tabLeg.tabQLeg

    # Compute Q_k(w)
    # Fill the initial values of Q_k(w)
    val_0_Q = log(1.0 - omg) - log(-1.0 - omg)
    val_1_Q = 2.0 + omg * val_0_Q
    # Compute the rest of the k values
    tabQLeg!(omg, val_0_Q, val_1_Q, tabQLeg)

    # Compute the value of H[-1,1,Re[w]]
    val_H = Heaviside(real(omg))

    # If we are outside of the interval [-1,1], we do not need to compute P_k(w) (it =0)
    if val_H == 0.0
        # Computing D_k(w)
        for k in 0:(struct_tabLeg.Ku-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] # Filling in tabDLeg
        end
    else
        # The container for P_k(w)
        tabPLeg = struct_tabLeg.tabPLeg
        # Initial values of P_k(w)
        val_0_P = 1.0 + 0.0im # Initial value for k = 0
        val_1_P = omg          # Initial value for k = 1
        tabPLeg!(omg, val_0_P, val_1_P, struct_tabLeg.Ku, tabPLeg) # Compute tabPLeg

        # Computing the D_k(w)
        for k in 0:(struct_tabLeg.Ku-1) # Loop over the Legendre indices
            tabDLeg[k+1] = tabQLeg[k+1] + 2.0im * π * val_H * tabPLeg[k+1] # Filling in tabDLeg
        end
    end
end
