module forcingFunctionsLog

export nend_forcing_func, nintermediate_forcing_func, n1_forcing_func


function nend_forcing_func(i, j, k, grid, clock, model_fields, indx)
    #Vᵢ = Main.find_Vi(Main.Rs[indx], Main.Hs[indx])
    Vᵢ = Main.Vs[indx]
    #uᵢ = find_steady_velocity(indx)
    # Compute derivative for nᵢ
    nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
    Fenc = Main.coll_freq_concentration_parameterisation(model_fields, indx)
    
    Fcoll = Fenc .* nᵢ

    # crystal size redistribution
    if Main.crystal_size_collision_redistribution == 1
        dn_coll = - Main.V₁/Vᵢ * Fcoll /Main.volume
    else
        dn_coll = Main.αVolconst[indx] * Fcoll /Main.volume
    end

    #if G₂[i, j, k] > 0 # growth
        #final_conc = G₂[i, j, k]*model_fields.n9[i, j, k]/(Vᵢ - V₂) + dn_coll[i, j, k] + udn_dz[i, j, k]
        #print(", udn/dz= ", udn_dz[1, 1, 1], ", ")
        #return @inbounds final_conc
    #else # melt
        #final_conc = (G₃[i, j, k]*model_fields.n10[i, j, k])/(Vᵢ - V₂) + dn_coll[i, j, k] + udn_dz[i, j, k]
        #return @inbounds  final_conc
    #end
    return dn_coll[i, j, k]

end

function nintermediate_forcing_func(i, j, k, grid, clock, model_fields, indx)
    #ρ = find_density(model_fields.T, model_fields.S)
    #uᵢ = find_steady_velocity(indx)

    # Compute derivative for nᵢ
    nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
    #nᵢ = max.(nᵢ, 0)

    # Extend velocity array (avoid index errors)
    #udn_dz = -nᵢ .* uᵢ /depth #use the faces below

    # Access nᵢ₋₁ and nᵢ₊₁ safely
    #n_im1 = getfield(model_fields, Symbol("n$(indx-1)"))  # Field nᵢ₋₁
    #n_ip1 = getfield(model_fields, Symbol("n$(indx+1)"))  # Field nᵢ₊₁

    # Compute growth rates and Vi values dynamically
    #Gᵢ = find_growth_rate(model_fields.T, Main.Rs[indx], Main.Hs[indx])
    #G_im1 = find_growth_rate(model_fields.T, Main.Rs[indx-1], Main.Hs[indx-1])
    #G_ip1 = find_growth_rate(model_fields.T, Main.Rs[indx+1], Main.Hs[indx+1])

    #Vᵢ = Main.find_Vi(Main.Rs[indx], Main.Hs[indx])
    Vᵢ = Main.Vs[indx]
    #V_im1 = find_Vi(Main.Rs[indx-1], Main.Hs[indx-1])
    #V_ip1 = find_Vi(Main.Rs[indx+1], Main.Hs[indx+1])

    # get collision frquency
    Fenc = Main.coll_freq_concentration_parameterisation(model_fields, indx)
    Fcoll = Fenc .* nᵢ

    # find Fcollᵦ
    if Main.crystal_size_collision_redistribution == 2
        βlist = Main.βs[indx]
        βVolconsts = Main.βVolconst[indx]
        Fcollβsum =  zeros(1, 1, Main.numz)
        for βindx in eachindex(βlist)
            β = βlist[βindx]
            nᵦ = getfield(model_fields, Symbol("n$β"))  # Dynamically get field `nᵢ`
            Fencᵦ = Main.coll_freq_concentration_parameterisation(model_fields, β)
            Fcollᵦ = Fencᵦ .* nᵦ
            Fcollβsum = Fcollβsum .+ βVolconsts[βindx] * Fcollᵦ
        end

        glist = Main.gs[indx]
        gVolconsts = Main.gVolconst[indx]
        Fcollgsum =  zeros(1, 1, Main.numz)
        for gindx in eachindex(glist)
            gpos = glist[gindx]
            nᵧ = getfield(model_fields, Symbol("n$gpos"))  # Dynamically get field `nᵢ`
            Fencᵧ = Main.coll_freq_concentration_parameterisation(model_fields, gpos)
            Fcollᵧ = Fencᵧ .* nᵧ
            Fcollgsum = Fcollgsum .+ gVolconsts[gindx] * Fcollᵧ
        end
    end



    # crystal size redistribution
    if Main.crystal_size_collision_redistribution == 1
        dn_coll = - Main.V₁/Vᵢ * Fcoll/Main.volume
    else
        dn_coll = (Main.αVolconst[indx] * Fcoll .+ Fcollβsum .+ Fcollgsum) / Main.volume
    end
    #print("n2", maximum(dn_coll))

    # Apply logic for growth and melting
    #if G_im1[i, j, k] > 0  # Growth case
    #    final_conc = - (Gᵢ[i, j, k] * nᵢ[i, j, k] / (V_ip1 - Vᵢ) - G_im1[i, j, k] * n_im1[i, j, k] / (Vᵢ - V_im1)) + dn_coll[i, j, k] + udn_dz[i, j, k]
    #    return @inbounds final_conc
    #else  # Melt case
    #    final_conc = - (G_ip1[i, j, k] * n_ip1[i, j, k] / (V_ip1 - Vᵢ) - Gᵢ[i, j, k] * nᵢ[i, j, k] / (Vᵢ - V_im1)) + dn_coll[i, j, k] + udn_dz[i, j, k]
    #    return @inbounds final_conc
    #end
    return @inbounds dn_coll[i, j, k]
end


function n1_forcing_func(i, j, k, grid, clock, model_fields, indx)
    #nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
    #nᵢ = max.(nᵢ, 0)
    #udn_dz = -nᵢ .* uᵢ /depth #use the faces below

    # growth term
    #G₁ = find_growth_rate(model_fields.T, Main.Rs[indx], Main.Hs[indx])
    #G₂ = find_growth_rate(model_fields.T, Main.Rs[indx+1], Main.Hs[indx+1])
    #Main.V₁ = find_Vi(Main.Rs[indx], Main.Hs[indx])
    #V₂ = find_Vi(Main.Rs[indx+1], Main.Hs[indx+1])

    # get collision frquency
    Fcollsum =  zeros(1, 1, Main.numz)
    for jindx in 2:lastindex(Main.Rs)
        Fenc = zeros(1, 1, Main.numz)
        nⱼ = getfield(model_fields, Symbol("n$jindx"))  # Dynamically get field `nᵢ`
        Fenc = Main.coll_freq_concentration_parameterisation(model_fields, jindx)
        Fcollsum = Fcollsum .+ Fenc .* nⱼ
    end  

    
    if Main.crystal_size_collision_redistribution == 1
        dn_coll = Fcollsum / Main.volume
    else
         # find Fcollᵦ
        if Main.crystal_size_collision_redistribution == 2
            βlist = Main.βs[indx]
            βVolconsts = Main.βVolconst[indx]
            Fcollβsum =  zeros(1, 1, Main.numz)
            for βindx in eachindex(βlist)
                β = βlist[βindx]
                nᵦ = getfield(model_fields, Symbol("n$β"))  # Dynamically get field `nᵢ`
                Fencᵦ = Main.coll_freq_concentration_parameterisation(model_fields, β)
                Fcollᵦ = Fencᵦ .* nᵦ
                Fcollβsum = Fcollβsum .+ βVolconsts[βindx] * Fcollᵦ
            end
        end
        dn_coll = (Main.ζ * Fcollsum .+ Fcollβsum) / Main.volume
    end

    #print("n1", maximum(dn_coll))

    #if G₁[i, j, k] > 0 # growth
    #    final_conc =   - G₁[i, j, k]*model_fields.n1[i, j, k]/(V₂ - Main.V₁) + dn_coll[i, j, k] + udn_dz[i, j, k]
    #    return @inbounds final_conc
    #else # melt
    #    final_conc = - (G₂[i, j, k]*model_fields.n2[i, j, k]/(V₂ - Main.V₁) - G₁[i, j, k]*model_fields.n1[i, j, k]/Main.V₁) + dn_coll[i, j, k] + udn_dz[i, j, k]
        
    return @inbounds  dn_coll[i, j, k]
    #end
end

end
