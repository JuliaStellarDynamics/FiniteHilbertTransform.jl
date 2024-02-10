# Adding the package src to the load path
push!(LOAD_PATH,"../src/")

using Documenter, FiniteHilbertTransform

makedocs(sitename = "FiniteHilbertTransform.jl",
         pages=[
                "Home" => "index.md",
                "Functions" => "functions.md",
                "Installation" => "installation.md",
                "Example" => "example.md"
               ],
         format = Documenter.HTML(prettyurls=false))

deploydocs(repo="github.com/JuliaStellarDynamics/FiniteHilbertTransform.jl",devbranch="documentation")