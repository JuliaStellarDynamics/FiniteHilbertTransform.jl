# FiniteHilbertTransform.jl

*Computation of Hilbert Transform with finite boundaries and Landau's prescription.*

While both Legendre and Chebyshev methods are implemented, we recommend using Legendre techniques for integration.

FiniteHilbertTransform.jl's core functionality precomputes the Hilbert-transformed Legendre functions \( Q_k(w) \) for a given complex frequency \( omg \). The Hilbert transform is defined as \( Q_k(w) = \int_{-1}^{1} \frac{P_k(u)}{u - w} du \), where \( P_k(u) \) is the Legendre function of the first kind. It is important to note that \( Q_k(w) = -2 q_k(w) \) for real values of \( w \), where \( q_k(w) \) represents the Legendre functions of the second kind.

