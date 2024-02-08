
# FiniteHilbertTransform.jl

[![image](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/FiniteHilbertTransform.jl/)

**FiniteHilbertTransform.jl** is a Julia package designed to compute the finite version of the Hilbert transformations. This toolbox is inspired by Tricomi's work on the finite Hilbert transform from 1957. In the context of gravitational dynamics, the finite Hilbert transform may be used as a scheme for analytic continuation to the lower half of the complex plane. See [Fouvry & Prunet (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2443F/abstract), or [Petersen et al. (2023)](https://ui.adsabs.harvard.edu/abs/2023arXiv231110630P/abstract) for details.

---
## Installation

**FiniteHilbertTransform** is currently unregistered[^1]. To add it to your julia[^2] registry, follow these steps:

1. **Add Package:** Use the package manager and execute the following command inside julia:
    ```julia
    add "git@github.com:JuliaStellarDynamics/FiniteHilbertTransform.git"
    ```
or at the command line

    ```bash
    $ julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl.git")'
    ```

2. **Verify Version:** Confirm the current version with `status FiniteHilbertTransform` in the julia package manager.

3. **Import Package:** Import the package in your julia environment with `import FiniteHilbertTransform`.

## Working from Source

Alternatively, work directly from the codebase:

1. **Activate Environment:** In the main directory of the package, enter the Julia environment using `julia`.

2. **Access Package Manager:** Inside the Julia environment, open the package manager with `]`.

3. **Activate Project:** Activate the project using `activate .`. For added safety, resolve dependencies using `resolve` to check for updates.

4. **Return to Julia Interpreter:** Exit the package manager with `[backspace]`. You are now equipped with the latest package version.

5. **Import Package:** Import the package by typing `using FiniteHilbertTransform` in the Julia interpreter.

Alternately[^3], you may clone the repository wherever you want and create a local environment (or project) by running:
```
$ git clone https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl.git
$ cd FiniteHilbertTransform.jl
$ julia --project=. -e 'using Pkg; Pkg.precompile()'
```

Note: If you are using a new Julia interpreter, you might need to download additional packages. Use the following command:
```julia
using(Pkg)
Pkg.instantiate()
```

---
## Quick use test

An introductory non-trivial example is given in `examples/run_plasma.jl`. This script will recreate Figure E1 from [Fouvry & Prunet (2021).](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2443F/abstract)

If you installed the library using the first (global) install option, just download this example [file](https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl/blob/main/examples/run_plasma.jl) from the github repository.

Run the code with the following command[^4]:
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


[^1]: For detailed instructions, check [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

[^2]:If you are new to `julia`, install the latest version by running this in your terminal: `$ curl -fsSL https://install.julialang.org | sh`. If you are on Windows or run into problems with `curl`-based installation, please visit [this website](https://julialang.org/downloads/).

[^3]: Note that if you use this install option you will always need to run codes in the project context by adding the option `--project=/path/to/FiniteHilbertTransform.jl` after `julia`. The library will not be accessible in your global julia context.

[^4]: Do not forget the option `--project=/path/to/FiniteHilbertTransform.jl` after `julia` if you installed the library locally.