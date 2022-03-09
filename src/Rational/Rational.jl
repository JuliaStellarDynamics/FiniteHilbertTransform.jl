# In this section, we perform the rational function
# approximation as in Weinberg

struct struct_Rational_type
    #=
    ##################################################
    # Defining a structure to perform the recurrences
    # for the rational approximations
    ##################################################
    =#
    taba    ::Array{Complex{Float64},1} # Static container for the values of the coefficients `a' in the continued fractions
    tabrho  ::Array{Complex{Float64},1} # Static container for the values of rho(w_1,..,w_i) (1<=i<=nbGrid_Rational)
    tabrho_1::Array{Complex{Float64},1} # Temporary container for rho(w_i,...,w_(i+k)) (i varies)
    tabrho_2::Array{Complex{Float64},1} # Temporary container for rho(w_i,...,w_(i+k)) (i varies)
    tabrho_3::Array{Complex{Float64},1} # Temporary container for rho(w_i,...,w_(i+k)) (i varies)
end

function struct_Rational_create(nbGrid)
    # Function that creates a struct_Rational
    return struct_Rational_type(zeros(Complex{Float64},nbGrid),
                                zeros(Complex{Float64},nbGrid),
                                zeros(Complex{Float64},nbGrid),
                                zeros(Complex{Float64},nbGrid),
                                zeros(Complex{Float64},nbGrid)) # Creating a struct_Rational with all the arrays set to 0.0+0.0im
end

function tabrho_Rational!(tabw::Array{Complex{Float64},1},
                          tabf::Array{Complex{Float64},1},
                          struct_Rational::struct_Rational_type)
    #=

    ##################################################
    # Function that for a given set of nodes (tabw)
    # and values (tabf) updates the array
    # of the reciprocal differences (tabrho)
    # @IMPROVE -- We could reduce the number of arrays reading
    ##################################################

    =#
    #####
    # Giving a name to the local variables
    tabrho   = struct_Rational.tabrho
    tabrho_1 = struct_Rational.tabrho_1
    tabrho_2 = struct_Rational.tabrho_2
    tabrho_3 = struct_Rational.tabrho_3
    #####
    nbGrid_Rational = size(tabrho)[1] # Number of elements in the grid
    #####
    # First, we fill in the value of rho(w_i)
    for i=1:nbGrid_Rational
        tabrho_1[i] = tabf[i] # Filling in the temporary array tabrho_1
    end
    tabrho[1] = tabrho_1[1] # Storing the value of rho(w_1)
    #####
    # Second, we fill in the initial value of rho(w_i,w_i+1)
    for i=1:(nbGrid_Rational-1) # ATTENTION, to the length of the loop
        tabrho_2[i] = (tabw[i] - tabw[i+1])/(tabf[i] - tabf[i+1])
    end
    tabrho[2] = tabrho_2[1] # Storing the value of rho(w_1,w_2)
    #####
    # Now, we can perform the recurrence
    for k=2:(nbGrid_Rational-1) # Number of times the recurrence has to be applied
        #####
        # Dealing with the generic term
        for i=1:(nbGrid_Rational-k)
            tabrho_3[i] = (tabw[i] - tabw[i+k])/(tabrho_2[i] - tabrho_2[i+1]) + tabrho_1[i+1] # Computing rho(w_i,w_{i+1},...,w_{i+k})
        end
        #####
        tabrho[k+1] = tabrho_3[1] # Storing the value of rho(w_1,w_2,...,w_{k})
        #####
        # Updating the temporary containers
        # for the next recurrence
        for i=1:nbGrid_Rational
            tabrho_1[i] = tabrho_2[i]
            tabrho_2[i] = tabrho_3[i]
        end
    end
end

function tab_Rational!(tabw::Array{Complex{Float64},1},
                       tabf::Array{Complex{Float64},1},
                       struct_Rational::struct_Rational_type)
    #=

     Function that computes the coefficients
     a_i (1<=i<=nbGrid_Rational)
     for a given rational approximation
     @IMPROVE -- we could slightly reduce the number of reads

    =#

    # Giving a name to the local variables
    tabrho = struct_Rational.tabrho
    taba   = struct_Rational.taba

    nbGrid_Rational = size(tabrho)[1] # Number of elements in the grid


    # First, we compute the values of rho(w_1,...,w_i) (1<=i<=nbGrid_Rational)
    tabrho_Rational!(tabw,tabf,struct_Rational)

    # Filling in the values of the two first coefficients
    taba[1] = tabrho[1]
    taba[2] = tabrho[2]


    # Generic expression for the coefficients
    for i=3:nbGrid_Rational
        taba[i] = tabrho[i] - tabrho[i-2]
    end
end


function evaluate_Rational(w::Complex{Float64},
                           tabw::Array{Complex{Float64},1},
                           taba::Array{Complex{Float64},1})
    #=

     Function to return the value of the continued
     fraction evaluated at a certain location
     ATTENTION, it is assumed that taba
     has already been updated

    =#

    nbGrid_Rational = size(taba)[1] # Number of elements in the grid

    f = taba[nbGrid_Rational] # Initialisation of the recurrence

    for i=(nbGrid_Rational-1):-1:1 # Loop over the coefficients. ATTENTION, this is a decreasing recurrence
        f = taba[i] + (w - tabw[i])/f # Update through the recurrence
    end

    return f # Output
end
