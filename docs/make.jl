# Adding the package src to the load path
push!(LOAD_PATH,"../src/")

using Documenter, FiniteHilbertTransform

makedocs(sitename = "FiniteHilbertTransform.jl",
         pages=[
                "Home" => "index.md",
                "Installation" => "installation.md",
                "Example" => "example.md",
                "Functions" => "functions.md"
               ],
         format = Documenter.HTML(prettyurls=false))

deploydocs(repo="github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl",devbranch="main")