"""
    tabLeg!_UNSTABLE(omg::Complex{Float64}, struct_tabLeg::LegendreFHT)

Computes the Legendre coefficients D_k(w) using the UNSTABLE integration method.

# Arguments
- `omg::Complex{Float64}`: Complex frequency.
- `struct_tabLeg::LegendreFHT`: LegendreFHT structure.

# Description
`tabLeg!_UNSTABLE` computes the Legendre coefficients D_k(w) using the UNSTABLE integration method. It computes the coefficients by first calculating the coefficients Q_k(w) and then setting D_k(w) equal to Q_k(w) for each k.

# Notes
- This function assumes that the container for D_k(w) (`struct_tabLeg.tabDLeg`) and Q_k(w) (`struct_tabLeg.tabQLeg`) are already initialized in the `struct_tabLeg`.

"""
function tabLeg!_UNSTABLE(omg::Complex{Float64}, struct_tabLeg::LegendreFHT)
    tabDLeg = struct_tabLeg.tabDLeg
    tabQLeg = struct_tabLeg.tabQLeg

    # Initial values of Q_k(w)
    val_0_Q = log(1.0 - omg) - log(-1.0 - omg) # Initial value for k = 0
    val_1_Q = 2.0 + omg * val_0_Q             # Initial value for k = 1

    tabQLeg!(omg, val_0_Q, val_1_Q, tabQLeg) # Computing the tabQLeg

    # Computing D_k(w)
    for k in 0:(struct_tabLeg.Ku-1)
        tabDLeg[k+1] = tabQLeg[k+1] # Filling in tabDLeg
    end
end
