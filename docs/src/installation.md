# Installation

**FiniteHilbertTransform** is currently unregistered. To add it to your julia registry, follow these steps:


1. **Read Documentation:** For detailed instructions, check [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

2. **Add Package:** Use the package manager and execute the following command:
    ```julia
    add "git@github.com:michael-petersen/FiniteHilbertTransform.git"
    ```
or at the command line

    ```bash
    $ julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl.git")'
    ```

3. **Handling Git Keys:** If you encounter Git key errors, register your private key using the Julia shell prompt (access with `;`), and point to your private key:
    ```julia
    ssh-add ~/.ssh/id_rsa
    ```


4. **Verify Version:** Confirm the current version with `status FiniteHilbertTransform` in the julia package manager.

5. **Import Package:** Import the package in your julia environment with `import FiniteHilbertTransform`.

## Working from Source

Alternatively, work directly from the codebase:

1. **Activate Environment:** In the main directory of the package, enter the Julia environment using `julia`.

2. **Access Package Manager:** Inside the Julia environment, open the package manager with `]`.

3. **Activate Project:** Activate the project using `activate .`. For added safety, resolve dependencies using `resolve` to check for updates.

4. **Return to Julia Interpreter:** Exit the package manager with `[backspace]`. You are now equipped with the latest package version.

5. **Import Package:** Import the package by typing `using FiniteHilbertTransform` in the Julia interpreter.

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


For detailed instructions, check [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages).

If you are new to `julia`, install the latest version by running this in your terminal: `$ curl -fsSL https://install.julialang.org | sh`. If you are on Windows or run into problems with `curl`-based installation, please visit [this website](https://julialang.org/downloads/).

Note that if you use this install option you will always need to run codes in the project context by adding the option `--project=/path/to/FiniteHilbertTransform.jl` after `julia`. The library will not be accessible in your global julia context.

Do not forget the option `--project=/path/to/FiniteHilbertTransform.jl` after `julia` if you installed the library locally.