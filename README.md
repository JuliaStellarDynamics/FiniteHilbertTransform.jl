
# FiniteHilbertTransform.jl

[![Test Builds](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/github/JuliaStellarDynamics/FiniteHilbertTransform.jl/graph/badge.svg?token=LO51PVUTU1)](https://codecov.io/github/JuliaStellarDynamics/FiniteHilbertTransform.jl)
[![image](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/FiniteHilbertTransform.jl/)

**FiniteHilbertTransform.jl** is a Julia package designed to compute the finite version of the Hilbert transformations. This toolbox is inspired by Tricomi's work on the finite Hilbert transform from 1957. In the context of gravitational dynamics, the finite Hilbert transform may be used as a scheme for analytic continuation to the lower half of the complex plane. See [Fouvry & Prunet (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2443F/abstract), or [Petersen et al. (2024)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract) for details.

---
## Installation

Install Julia by following the instructions at [julialang.org/downloads/](https://julialang.org/downloads/).

To invoke Julia in the Terminal, you need to make sure that the `julia` command-line program is in your `PATH`. 
See [here](https://julialang.org/downloads/platform/#optional_add_julia_to_path) for detailed instructions.

Once Julia installed, obtain the `FiniteHilbertTransform.jl` library[^1][^2] and compile it by running:
    ```
    julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl.git")'
    ```

---
## Quickstart

An introductory example is given in `examples/run_plasma.jl`. This script will recreate Figure E1 from [Fouvry & Prunet (2021).](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2443F/abstract)

Download the file by running:
```
wget https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/blob/main/examples/run_plasma.jl
```
Run the code with the following command:
```
$ julia /path/to/run_plasma.jl
```

This example will first install some required libraries (`Plots`, `ArgParse`) as needed. These installations might take a few minutes when first called.

The resulting plot will be created in the same folder where you executed the script, and will be called `plasmademo.png`.

![`Plasma Demonstration`](examples/plasmademo.png)

---
### Interactive notebook

If you prefer interactive Jupyter notebooks, you will need to install `IJulia` following these [instructions](https://github.com/JuliaLang/IJulia.jl).

The interactive introduction example is then given in `examples/run_plasma.ipynb`.


---
### Without installing Julia

For those who prefer not to install julia locally, we also provide a Google colab version[^3] that may be run in the cloud. [See here](https://colab.research.google.com/drive/1p4lX5ot5-kKSnIo1XLFchsiOWGQUxEhR).


---
## Documentation and usage

To get more familiar with the content of the library and start and design your own use case, you may want to visit the [documentation](https://juliastellardynamics.github.io/FiniteHilbertTransform.jl/).


----------------------------

## Authors

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk

Mathieu Roule -  @MathieuRoule - roule@iap.fr



[^1]: The library is also easy to uninstall: remove the package from the environment by running
```
julia -e 'using Pkg; Pkg.rm("FiniteHilbertTransform");'
```

[^2]: By default, packages are added to the default environment at ~/.julia/environments/v1.#. It is however easy to create other, independent, projects. If you want to install the `FiniteHilbertTransform` package in a different/test environment, first create a folder to host the environment files (Project.toml and Manifest.toml which will be created later on). Then, for every command line invoking Julia, use `julia --project=/path/to/my_env` instead of `julia` alone.

[^3]: This notebook is not maintained as a priority. We would recommand you install Julia on your machine to test the library locally.
