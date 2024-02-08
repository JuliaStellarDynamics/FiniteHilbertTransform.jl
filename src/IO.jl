using HDF5 # To have access to .hf5 formats

function dump_tabIminusXi(namefile::String,
                          tabomega::Vector{ComplexF64},
                          tabIminusXi::Vector{ComplexF64})
    #=dump_tabIminusXi

    dump the primary mode table

    =#
    file = h5open(namefile,"w") # Opening the file
    write(file,"tabomega_real",real.(tabomega)) # Writing the REAL part of the complex frequencies
    write(file,"tabomega_imag",imag.(tabomega)) # Writing the IMAG part of the complex frequencies
    write(file,"tabIminusXi_real",real.(tabIminusXi)) # Writing the REAL part of det[I-Xi]
    write(file,"tabIminusXi_imag",imag.(tabIminusXi)) # Writing the IMAG part of det[I-Xi]
    close(file) # Closing the file
end
