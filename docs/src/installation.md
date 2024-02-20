
# Installation

Install Julia by following the instructions at [julialang.org/downloads/](https://julialang.org/downloads/).

To invoke Julia in the Terminal, you need to make sure that the `julia` command-line program is in your `PATH`. 
See [here](https://julialang.org/downloads/platform/#optional_add_julia_to_path) for detailed instructions.

Once Julia installed, obtain the `FiniteHilbertTransform.jl` library[^1][^2] and compile it by running:
    ```
    julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl.git")'
    ```

---
## From source

Alternately, you may clone the repository wherever you want and create a local environment (or project) by running:
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


[^1]: The library is also easy to uninstall: remove the package from the environment by running
```
julia -e 'using Pkg; Pkg.rm("FiniteHilbertTransform");'
```

[^2]: By default, packages are added to the default environment at ~/.julia/environments/v1.#. It is however easy to create other, independent, projects. If you want to install the `FiniteHilbertTransform` package in a different/test environment, first create a folder to host the environment files (Project.toml and Manifest.toml which will be created later on). Then, for every command line invoking Julia, use `julia --project=/path/to/my_env` instead of `julia` alone.
