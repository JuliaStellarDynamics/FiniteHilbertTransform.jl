using HDF5 # To have access to .hf5 formats

"""
    dump_tabIminusXi(namefile::String, tabomega::Vector{Complex{Float64}}, tabIminusXi::Vector{Complex{Float64}})

Dump the primary mode table to an HDF5 file.

# Arguments
- `namefile::String`: Name of the HDF5 file to write.
- `tabomega::Vector{Complex{Float64}}`: Vector of complex frequency values.
- `tabIminusXi::Vector{Complex{Float64}}`: Vector of complex values representing det[I - Xi].

# Description
`dump_tabIminusXi` writes the complex frequency values and det[I - Xi] values to an HDF5 file specified by `namefile`.

# Example
```julia
tabomega = [1.0 + 2.0im, 2.0 + 3.0im, 3.0 + 4.0im]
tabIminusXi = [0.5 + 0.1im, 0.6 + 0.2im, 0.7 + 0.3im]
dump_tabIminusXi("data.h5", tabomega, tabIminusXi)
```

#Notes
- The length of tabomega and tabIminusXi must be the same.
- The resulting HDF5 file will contain datasets named "tabomega_real", "tabomega_imag", "tabIminusXi_real", and "tabIminusXi_imag", storing the real and imaginary parts of tabomega and tabIminusXi, respectively.


"""
function dump_tabIminusXi(namefile::String,
                          tabomega::Vector{ComplexF64},
                          tabIminusXi::Vector{ComplexF64})

    file = h5open(namefile,"w") # Opening the file
    write(file,"tabomega_real",real.(tabomega)) # Writing the REAL part of the complex frequencies
    write(file,"tabomega_imag",imag.(tabomega)) # Writing the IMAG part of the complex frequencies
    write(file,"tabIminusXi_real",real.(tabIminusXi)) # Writing the REAL part of det[I-Xi]
    write(file,"tabIminusXi_imag",imag.(tabIminusXi)) # Writing the IMAG part of det[I-Xi]
    close(file) # Closing the file
end
