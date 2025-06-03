
module collisionFrequencyFunctions

export coll_freq_concentration_parameterisation_num_2, coll_freq_concentration_parameterisation_num_1

function get_Fefficiency_ndensity(indx, Rindx)
    #collision_efficiency = Main.collision_efficiency
    Feff = Main.collision_efficiency[indx, Rindx]
    return Feff
end

function get_Fefficiency_nsum(indx)
    #collision_efficiency = Main.collision_efficiency
    Feff = Main.collision_efficiency[indx, indx]
    return Feff
end


function coll_freq_concentration_parameterisation_num_2(model_fields, indx)
    Fcoll = zeros(1, 1, Main.numz)
    for Rindx in eachindex(Main.Rs)
        nRindx = getfield(model_fields, Symbol("n$Rindx"))
        Fcoll .+= Main.encounterFrequencyConstants_matrix[indx, Rindx]*nRindx
        #Fcoll .+= Main.get_Fenc(indx, Rindx, nRindx) #* get_Fefficiency_ndensity(indx, Rindx)

    end   
    return Fcoll
end

function coll_freq_concentration_parameterisation_num_1(model_fields, indx)
    # get collision frquency
    Fcoll = zeros(1, 1, Main.numz)
    nTotal = zeros(1, 1, Main.numz)
    for Rindx in eachindex(Main.Rs)
        nRindx = getfield(model_fields, Symbol("n$Rindx"))
        nTotal .+= nRindx
    end
    ntot = min.(nTotal, Main.nmax)
    Fcoll = Main.encounterFrequencyConstants_matrix[indx]*ntot
    #Fcoll = Main.get_Fenc_ndensity(indx, ntot) #* get_Fefficiency_nsum(indx)

    return Fcoll
end

end