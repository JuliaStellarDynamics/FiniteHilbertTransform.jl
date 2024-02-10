## Quick use test

An introductory non-trivial example is given in `examples/run_plasma.jl`. This script will recreate Figure E1 from [Fouvry & Prunet (2021).](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2443F/abstract)

If you installed the library using the first (global) install option, just download this example [file](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/blob/main/examples/run_plasma.jl) from the github repository.

Run the code with the following command:
```
$ julia /path/to/run_plasma.jl
```

This example will first install some required libraries (`Plots`, `ArgParse`). These installations might take a few minutes when first called.

The resulting plot will be created in the same folder as the test code under the name `plasmademo.png`.

![`Plasma Demonstration`](../../examples/plasmademo.png)
