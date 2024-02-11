
# FiniteHilbertTransform.jl

[![Test Builds](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/actions/workflows/benchmark.yml/badge.svg)](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/actions/workflows/benchmark.yml)
[![image](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/FiniteHilbertTransform.jl/)

**FiniteHilbertTransform.jl** is a Julia package designed to compute the finite version of the Hilbert transformations. This toolbox is inspired by Tricomi's work on the finite Hilbert transform from 1957. In the context of gravitational dynamics, the finite Hilbert transform may be used as a scheme for analytic continuation to the lower half of the complex plane. See [Fouvry & Prunet (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2443F/abstract), or [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract) for details.

---
## Installation

Install Julia by following the instructions at [julialang.org/downloads/](https://julialang.org/downloads/).

To invoke Julia in the Terminal, you need to make sure that the `julia` command-line program is in your `PATH`. 
See [here](https://julialang.org/downloads/platform/#optional_add_julia_to_path) for detailed instructions.

Once Julia installed, clone the `FiniteHilbertTransform.jl` library and precompile it by running:
```
git clone https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl.git
cd FiniteHilbertTransform.jl
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

---
## Quick use test

An introductory non-trivial example is given in `examples/run_plasma.jl`. This script will recreate Figure E1 from [Fouvry & Prunet (2021).](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2443F/abstract)

If you installed the library using the first (global) install option, just download this example [file](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/blob/main/examples/run_plasma.jl) from the github repository.

Run the code with the following command[^1]:
```
$ julia /path/to/run_plasma.jl
```

This example will first install some required libraries (`Plots`, `ArgParse`). These installations might take a few minutes when first called.

The resulting plot will be created in the same folder as the test code under the name `plasmademo.png`.

![`Plasma Demonstration`](examples/plasmademo.png)

### Interactive notebooks

If you prefer interactive Jupyter notebooks, you will need to install `IJulia` following these [instructions](https://github.com/JuliaLang/IJulia.jl).

The interactive introduction example is then given in `examples/run_plasma.ipynb`.

For those who prefer not to install julia locally, we also provide a Google colab version that may be run in the cloud. [See here](https://colab.research.google.com/drive/1p4lX5ot5-kKSnIo1XLFchsiOWGQUxEhR).

---
## Documentation and usage

To get more familiar with the content of the library and start and design your own use case, you may want to visit the [documentation](https://juliastellardynamics.github.io/FiniteHilbertTransform.jl/).


-----------------------------

### Authors

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk

Mathieu Roule -  @MathieuRoule - roule@iap.fr



[^1]: Do not forget the option `--project=/path/to/FiniteHilbertTransform.jl` after `julia` if you installed the library locally.