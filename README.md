
# PerturbPlasma.jl

`PerturbPlasma.jl` is a package written in Julia as a playground for plasma dispersion calculations.

-----------------------------

## Quick activate

In the main directory where you the package lives, enter the Julia environment (`julia`), then the package manager (`]`), then activate (`activate .`). To be extra safe, you can `resolve` to check for updates. Then return to the Julia interpreter (`[backspace]`): you are good to go with the latest version of the package! Import by typing `using PerturbPlasma` into the Julia interpreter.

-----------------------------

## Testing Plasma Techniques

`0` will do the setup work of computing the coefficients at K_u Legendre nodes for a plasma model defined by qself and xmax.

`compute_tabIminusXi(tabomega,taba,xmax,struct_tabLeg)` will perform the Legendre integration (struct_tabLeg) over the coefficients (taba) at an array of complex frequencies (tabomega) for a model parameterised by frequency xmax.

-----------------------------

## Author

Mike Petersen -  @michael-petersen - petersen@iap.fr
