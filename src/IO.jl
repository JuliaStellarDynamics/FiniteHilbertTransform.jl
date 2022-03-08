using HDF5 # To have access to .hf5 formats

function dump(namefile::String,
              tabomega::Vector{Complex{Float64}},
              tabIminusXi::Vector{Complex{Float64}})
    file = h5open(namefile,"w") # Opening the file
    #write(file,"tabOmega",tabOmega) # Writing tabOmega to the file
    #write(file,"tabEta",tabEta) # Writing tabEta to the file
    write(file,"tabomega_real",real.(tabomega)) # Writing the REAL part of the complex frequencies
    write(file,"tabomega_imag",imag.(tabomega)) # Writing the IMAG part of the complex frequencies
    write(file,"tabIminusXi_real",real.(tabIminusXi)) # Writing the REAL part of det[I-Xi]
    write(file,"tabIminusXi_imag",imag.(tabIminusXi)) # Writing the IMAG part of det[I-Xi]
    close(file) # Closing the file
end
