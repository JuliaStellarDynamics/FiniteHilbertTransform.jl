
# PerturbPlasma.jl

`PerturbPlasma.jl` is a package written in Julia as a playground for plasma dispersion calculations.

-----------------------------

### Quick activate

`PerturbPlasma` is (currently) unregistered, and as such if you would like to add it to your Julia registry, read [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages). Short version: when in the package manager, `add "git@github.com:michael-petersen/PerturbPlasma.git"`. (If you are getting an error about git keys, you will need to register your private key using the julia shell prompt (access with `;`), and then pointing at your private key: `ssh-add ~/.ssh/id_rsa`.) You can verify the current version with `status PerturbPlasma` in the package manager. Then, `import PerturbPlasma` and you are ready to go.

Another option, to work directly from a codebase: in the main directory where you the package lives, enter the Julia environment (`julia`), then the package manager (`]`), then activate (`activate .`). To be extra safe, you can `resolve` to check for updates. Then return to the Julia interpreter (`[backspace]`): you are good to go with the latest version of the package! Import by typing `using PerturbPlasma` into the Julia interpreter. You may also need to download some packages if you are using a new Julia interpreter: try `using(Pkg);Pkg.instantiate()`.

-----------------------------

### Testing Plasma Techniques (Legendre)

`tabuGLquad,tabwGLquad = PerturbPlasma.tabuwGLquad(K_u)` will call out for the u positions (tabuGLquad) and corresponding weights (tabwGLquad) at K_u points.

`setup_legendre_integration(K_u,qself,xmax,parallel)` will do the setup work of computing the coefficients at K_u Legendre nodes for a plasma model defined by qself and xmax.

`compute_tabIminusXi(tabomega,taba,xmax,struct_tabLeg)` will perform the Legendre integration (struct_tabLeg) over the coefficients (taba) at an array of complex frequencies (tabomega) for a model parameterised by frequency xmax.

`get_Legendre_IminusXi(omega,taba,xmax,struct_tabLeg)` will perform the Legendre integration (struct_tabLeg) over the coefficients (taba) for a single complex frequency (omega) for a model parameterised by frequency xmax.

With these commands, you can perform a manual search for the zeros of I - Xi, the response matrix. The location of zeros indicates the presence of a mode.
```using PerturbPlasma
tabaL,structL = setup_legendre_integration(75,0.5,20.,false)
test_ninepointsL(tabaL,20.,structL)
```

-----------------------------

### Testing Plasma Techniques (Chebyshev)

`setup_chebyshev_integration(K_u,qself,xmax,parallel)` will do the setup work of computing the coefficients at K_u Legendre nodes for a plasma model defined by qself and xmax.

`compute_tabIminusXi(tabomega,taba,xmax,method)` will perform the Chebyshev integration (using specified method) over the coefficients (taba) at an array of complex frequencies (tabomega) for a model parameterised by frequency xmax.

`get_Chebyshev_IminusXi(omega,taba,xmax)` will perform the Chebyshev integration over the coefficients (taba) for a single complex frequency (omega) for a model parameterised by frequency xmax.

Complete demo:
```using PerturbPlasma
tabaC = setup_chebyshev_integration(75,0.5,20.,false)
test_ninepointsC(tabaC,20.)
```

-----------------------------

### Author

Mike Petersen -  @michael-petersen - petersen@iap.fr
