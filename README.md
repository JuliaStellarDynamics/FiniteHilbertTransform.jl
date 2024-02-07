
# FiniteHilbertTransform.jl
## Version 0.9

[![image](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/FiniteHilbertTransform.jl/)

**FiniteHilbertTransform.jl** is a Julia package designed to compute finite Hilbert transformations. This toolbox is inspired by Tricomi's work from 1957, offering powerful capabilities despite its vintage influence.

## Quick Activation

**FiniteHilbertTransform** is currently unregistered. To add it to your Julia registry, follow these steps:

1. **Read Documentation:** For detailed instructions, check [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

2. **Add Package:** Use the package manager and execute the following command:
    ```julia
    add "git@github.com:michael-petersen/FiniteHilbertTransform.git"
    ```

3. **Handling Git Keys:** If you encounter Git key errors, register your private key using the Julia shell prompt (access with `;`), and point to your private key:
    ```julia
    ssh-add ~/.ssh/id_rsa
    ```

4. **Verify Version:** Confirm the current version with `status FiniteHilbertTransform` in the package manager.

5. **Import Package:** Import the package in your Julia environment with `import FiniteHilbertTransform`.

## Working from Source

Alternatively, work directly from the codebase:

1. **Activate Environment:** In the main directory of the package, enter the Julia environment using `julia`.

2. **Access Package Manager:** Inside the Julia environment, open the package manager with `]`.

3. **Activate Project:** Activate the project using `activate .`. For added safety, resolve dependencies using `resolve` to check for updates.

4. **Return to Julia Interpreter:** Exit the package manager with `[backspace]`. You are now equipped with the latest package version.

5. **Import Package:** Import the package by typing `using FiniteHilbertTransform` in the Julia interpreter.

Note: If you are using a new Julia interpreter, you might need to download additional packages. Use the following command:
```julia
using(Pkg)
Pkg.instantiate()
```

**FiniteHilbertTransform.jl** provides efficient and reliable tools for finite Hilbert transformations, ensuring seamless integration with your Julia workflows.


This function precomputes the Hilbert-transformed Legendre functions \( Q_k(w) \) for a given complex frequency \( omg \). The Hilbert transform is defined as \( Q_k(w) = \int_{-1}^{1} \frac{P_k(u)}{u - w} du \), where \( P_k(u) \) is the Legendre function of the first kind. It is important to note that \( Q_k(w) = -2 q_k(w) \) for real values of \( w \), where \( q_k(w) \) represents the Legendre functions of the second kind.

### Testing Plasma Techniques (Legendre)

`tabuGLquad,tabwGLquad = FiniteHilbertTransform.tabuwGLquad(K_u)` will call out for the u positions (tabuGLquad) and corresponding weights (tabwGLquad) at K_u points.

`setup_legendre_integration(K_u,qself,xmax,parallel)` will do the setup work of computing the coefficients at K_u Legendre nodes for a plasma model defined by qself and xmax.

`compute_tabIminusXi(tabomega,taba,xmax,struct_tabLeg)` will perform the Legendre integration (struct_tabLeg) over the coefficients (taba) at an array of complex frequencies (tabomega) for a model parameterised by frequency xmax.

`get_Legendre_IminusXi(omega,taba,xmax,struct_tabLeg)` will perform the Legendre integration (struct_tabLeg) over the coefficients (taba) for a single complex frequency (omega) for a model parameterised by frequency xmax.

With these commands, you can perform a manual search for the zeros of I - Xi, the response matrix. The location of zeros indicates the presence of a mode.
```using FiniteHilbertTransform
tabaL,structL = setup_legendre_integration(75,0.5,20.,false)
test_ninepointsL(tabaL,20.,structL)
```

-----------------------------

### Testing Plasma Techniques (Chebyshev)

`setup_chebyshev_integration(K_u,qself,xmax,parallel)` will do the setup work of computing the coefficients at K_u Legendre nodes for a plasma model defined by qself and xmax.

`compute_tabIminusXi(tabomega,taba,xmax,method)` will perform the Chebyshev integration (using specified method) over the coefficients (taba) at an array of complex frequencies (tabomega) for a model parameterised by frequency xmax.

`get_Chebyshev_IminusXi(omega,taba,xmax)` will perform the Chebyshev integration over the coefficients (taba) for a single complex frequency (omega) for a model parameterised by frequency xmax.

Complete demo:
```using FiniteHilbertTransform
tabaC = setup_chebyshev_integration(75,0.5,20.,false)
test_ninepointsC(tabaC,20.)
```

-----------------------------

### Authors

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk

Mathieu Roule -  @MathieuRoule - roule@iap.fr
