"""
    tabLeg!_NEUTRAL(omg::Complex{Float64}, struct_tabLeg::LegendreFHT)

Fill in all the Legendre arrays for a NEUTRAL mode, i.e., Im[w] = 0.0.

# Arguments
- `omg::Complex{Float64}`: Complex frequency.
- `struct_tabLeg::LegendreFHT`: LegendreFHT structure.

# Description
`tabLeg!_NEUTRAL` fills in all the Legendre arrays for a NEUTRAL mode, where Im[w] = 0.0. It computes the Legendre coefficients D_k(w) and, if needed, Legendre polynomials P_k(w) for the given complex frequency.

"""
function tabLeg!_NEUTRAL(omg::Complex{Float64}, struct_tabLeg::LegendreFHT)
    # Simplify naming convention for the container for D_k(w) and Q_k(w)
    tabDLeg = struct_tabLeg.tabDLeg
    tabQLeg = struct_tabLeg.tabQLeg
    tabPLeg = struct_tabLeg.tabPLeg

    # Keep only the real part of omg
    romg = real(omg)

    # Initial values of Q_k(w):
    if romg == 1.0
        # Unique case, requires special seed values
        val_0_Q = 0.0 + 0.0im
        val_1_Q = 2.0 + 0.0im
    else
        # Initial value for k = 0
        val_0_Q = log(abs(1.0 - romg)) - log(abs(-1.0 - romg)) + 0.0im
        # Initial value for k = 1
        val_1_Q = 2.0 + romg * val_0_Q + 0.0im
    end

    # Compute tabQLeg
    tabQLeg!(omg, val_0_Q, val_1_Q, tabQLeg)

    # Compute P_k(w) if needed
    val_H = Heaviside(romg)

    # We are outside of the interval [-1,1], no need to compute P_k(w)
    if val_H == 0.0
        # Compute D_k(w)
        for k in 1:struct_tabLeg.Ku
            tabDLeg[k] = tabQLeg[k]
        end
    else
        # Initial values of P_k(w)
        val_0_P = 1.0 + 0.0im # Initial value for k = 0
        val_1_P = romg + 0.0im # Initial value for k = 1

        tabPLeg!(omg, val_0_P, val_1_P, struct_tabLeg.Ku, tabPLeg) # Compute tabPLeg

        # Compute D_k(w)
        for k in 1:struct_tabLeg.Ku
            tabDLeg[k] = tabQLeg[k] + im * π * val_H * tabPLeg[k]
        end
    end
end
