module SetupFunctions

export find_mean_collision_radius, find_collision_efficiency, find_β_α_Vol_constants, find_steady_velocity_via_iteration, find_Vi, find_β_α_Vol_constants_log, find_vgosink, find_vstokes, velocity_transition

# setup functions i.e. the areas and velocities
#define a matrix which gives the effective collision area when allowing for the different orientations, call it collision radius
function find_mean_collision_radius(Rs, aspect_ratio)
    effective_areas = zeros(length(Rs), length(Rs))
    for Riindx in eachindex(Rs)
        rᵢ = Rs[Riindx]
        for Rjindx in eachindex(Rs)
            rⱼ = Rs[Rjindx]
            effective_areas[Riindx, Rjindx] = 2*(rᵢ + rⱼ)*(1 + aspect_ratio)/(aspect_ratio * π)
        end
    end
    return effective_areas
end

function find_collision_efficiency(Rs)
    collision_efficiencies = [0, 0.1776, 0.2103, 0.4725, 0.5859, 0.5859]
    ratios = [1, 2, 5, 10, 100, 200]
    interp_linear = LinearInterpolation(ratios, collision_efficiencies)

    collision_efficiency = zeros(length(Rs), length(Rs))
    for Riindx in eachindex(Rs)
        rᵢ = Rs[Riindx]
        for Rjindx in eachindex(Rs)
            rⱼ = Rs[Rjindx]
            if rⱼ > rᵢ
                ratio = rⱼ/rᵢ
            else
                ratio = rᵢ/rⱼ
            end
            collision_efficiency[Riindx, Rjindx] = interp_linear(ratio)
            #
            
            #if ratio > 10
            #    collision_efficiency[Riindx, Rjindx] = 0.5859
            #elseif ratio > 5
            #    collision_efficiency[Riindx, Rjindx] = 0.4725
            #elseif ratio > 2
            #    collision_efficiency[Riindx, Rjindx] = 0.2103
            #else
            #    collision_efficiency[Riindx, Rjindx] = 0.1776
            #end
        end
    end

    return collision_efficiency
end

function find_β_α_Vol_constants(Vs, Vrem, V₁, ζ)
    # fint the values of size class the fractures into/into which a crystal fractures and the corresponding volume distribution

    # matrix of which crystal size the crystal breaks off into
    αs = zeros(length(Vs))
    for Vindx in eachindex(Vs)
        Vpostcoll = Vs[Vindx] - Vrem
        if Vpostcoll >= V₁
            αs[Vindx] = findlast(Vpostcoll .- Vs .>= 0) # find the biggest crystal smaller than the volume of the crystal that is left
        end
    end
    αs = Int.(αs)

    # matrix of which crystal breaking introduces crystals of this size
    βs = zeros(length(Vs))
    for Vindx in eachindex(Vs)
        try
            βs[Vindx] = findfirst(x -> x == Vindx, αs)
        catch
        end
    end
    βs = Int.(βs)
    # coefficient for crystal of size class i breaking
    αVolconst = zeros(length(Vs))
    for Vindx in eachindex(αVolconst)
        α = αs[Vindx]
        if α != 0
            αVolconst[Vindx] = -ζ * V₁/Vs[Vindx] * (1 + Vs[α]/(Vs[Vindx] - Vs[α]))
        end
    end
    # coefficient for crystal breaking into size class i
    βVolconst = zeros(length(Vs))
    for Vindx in eachindex(βVolconst)
        β = βs[Vindx]
        if β != 0
            βVolconst[Vindx] = ζ * V₁/(Vs[β] - Vs[Vindx])
        end
    end

    return αs, βs, αVolconst, βVolconst
end

function find_β_α_Vol_constants_log(Vs, Vrem, V₁, ζ)
    # fint the values of size class the fractures into/into which a crystal fractures and the corresponding volume distribution

    # matrix of which crystal size the crystal breaks off into (smallest size class)
    αs = zeros(length(Vs))
    for Vindx in eachindex(Vs)
        Vpostcoll = Vs[Vindx] - Vrem
        if Vpostcoll > V₁
            αs[Vindx] = findlast(Vpostcoll .- Vs .>= 0)
        else
            αs[Vindx] = 1
        end
    end
    αs = Int.(αs)

    # matrix of the larger size class in the sandwich
    ls = zeros(length(Vs))
    for Vindx in eachindex(Vs)
        Vpostcoll = Vs[Vindx] - Vrem
        if Vpostcoll > V₁
            ls[Vindx] = findfirst(Vs .- Vpostcoll .>= 0)
        else
            ls[Vindx] = Vindx # otherwise the other size class will be the same as before
        end
        if ls[Vindx] == Vindx
            ls[Vindx] = 0 # set it to zero if the larger size class is size class i
        end
    end
    ls = Int.(ls)

    # find the volume coefficient αVolconst for size class i, sVolconst for size class α and lVolconst for size class β
    αVolconst = zeros(length(Vs))
    sVolconst = zeros(length(Vs))
    lVolconst = zeros(length(Vs))
    for Vindx in eachindex(αVolconst)
        α = αs[Vindx]
        if α != 1
            l = ls[Vindx]
            if l == 0 # if the larger size class in the sandwich is of size class Vindx
                αVolconst[Vindx] = -ζ * V₁/(Vs[Vindx] - Vs[α])
                sVolconst[Vindx] = ζ * V₁/(Vs[Vindx] - Vs[α])
            else
                αVolconst[Vindx] = -1
                sVolconst[Vindx] = (Vs[l] - Vs[Vindx] + ζ * V₁)/(Vs[l] - Vs[α])
                lVolconst[Vindx] = (Vs[Vindx] - ζ * V₁ - Vs[α])/(Vs[l] - Vs[α])
            end
        else # just breaks into the smallest size class
            αVolconst[Vindx] = -1
            sVolconst[Vindx] = Vs[Vindx]/V₁ - ζ
        end
    end

    # when a crystal of size class βs[i] breaks, the smaller crystal is of size class i.
    #βVolconstlist is the volume introduced into size class i
    βs = Dict{Int, Vector{Int}}()
    βVolconst = Dict{Int, Vector{Float64}}()
    for Vindx in eachindex(Vs)
        alphs_list = findall(αs .== Vindx)
        βs[Vindx] = alphs_list
        println(", alphs_list, ", alphs_list)
        βVolconstlist = zeros(length(alphs_list))
        # coefficient for crystal breaking into size class i
        for aindx in eachindex(alphs_list)
            βVolconstlist[aindx] =  sVolconst[alphs_list[aindx]]
            
            #if Vs[Vindx] - Vrem < V₁
            #    βVolconstlist[aindx] = Vs[alphs_list[aindx]]/V₁ - ζ
            #    println("a, ", βVolconstlist[aindx])
            #else
            #    βVolconstlist[aindx] =  ζ * V₁/(Vs[alphs_list[aindx]] - Vs[Vindx])
            #    
            #end
        end
        βVolconst[Vindx] = βVolconstlist
    end

    gs = Dict{Int, Vector{Int}}()
    gVolconst = Dict{Int, Vector{Float64}}()
    for Vindx in eachindex(Vs)
        gs_list = findall(ls .== Vindx)
        gs[Vindx] = gs_list
        println(", gs_list, ", gs_list)
        gVolconstlist = zeros(length(gs_list))
        # coefficient for crystal breaking into size class i
        for gindx in eachindex(gs_list)
            gVolconstlist[gindx] =  lVolconst[gs_list[gindx]]
        end
        gVolconst[Vindx] = gVolconstlist
    end    

    return αs, βs, αVolconst, βVolconst, gs, gVolconst
end

function find_steady_velocity_via_iteration_old(Rs, ν, ρₒ, ρᵢ, grav, aspect_ratio)
    vginitial = range(start = 1e-8, stop = 1e-1, step = 1e-9)
    vgs = zeros(length(Rs))
    for Rindx in eachindex(Rs)
        rᵢ = Rs[Rindx]
        Reₚ = 2 * rᵢ * vginitial/ν
        logCD = 1.386 .- 0.892*log10.(Reₚ) .+ 0.111 * (log10.(Reₚ)).^2
        #CD = exp.(logCD)
        CD = 10 .^ (logCD)
        #print(vginitial.^2)
        #print(4*(ρₒ - ρᵢ)/ρₒ*grav * rᵢ/aspect_ratio * 1 ./CD)
        #print(maximum(4*(ρₒ - ρᵢ)/ρₒ * grav * rᵢ/aspect_ratio * 1 ./CD))
        #print(abs.(vginitial.^2 .- 4*(ρₒ - ρᵢ)/ρₒ*grav * rᵢ/aspect_ratio * 1 ./CD))
        vgindx = argmin(abs.(vginitial.^2 .- 4*(ρₒ - ρᵢ)/ρₒ*grav * rᵢ/aspect_ratio * 1 ./CD))
        vgs[Rindx] = vginitial[vgindx]
    end
    return vgs
end

function find_vgosink(Rs, ν, ρₒ, ρᵢ, grav, aspect_ratio)
    #vginitial = range(start = 1e-8, stop = 1e-1, step = 1e-9)
    vgs = zeros(length(Rs))
    for Rindx in eachindex(Rs)
        rᵢ = Rs[Rindx]
        if rᵢ < 6e-5
            vginitial = 10 .^ (range(start=log10(1e-16), stop=log10(1e-6), length=1000000))
        elseif rᵢ < 1e-4
            vginitial = 10 .^ (range(start=log10(1e-8), stop=log10(1e-5), length=10000000))
        else
            vginitial = 10 .^ (range(start=log10(1e-5), stop=log10(1), length=1000000))
        end
        Reₚ = 2 * rᵢ * vginitial/ν
        logCD = 1.386 .- 0.892*log10.(Reₚ) .+ 0.111 * (log10.(Reₚ)).^2
        #CD = exp.(logCD)
        CD = 10 .^ (logCD)
        vgindx = argmin(abs.(vginitial.^2 .- 4*(ρₒ - ρᵢ)/ρₒ*grav * rᵢ/aspect_ratio * 1 ./CD))
        vgs[Rindx] = vginitial[vgindx]
    end
    return vgs
end
        
function find_vstokes(Rs, ν, ρₒ, ρᵢ, grav, aspect_ratio)
    vstokes = (ρₒ - ρᵢ)/ρₒ * grav * π *(2*Rs) .^ 2 / (32 * aspect_ratio * ν)
    return vstokes
end

function velocity_transition(vRHS, vLHS, Re_gos)
    vmixed = zeros(length(vRHS))
    #pmid = argmin(abs.(Re_gos .- 0.1)) # 29 is where Re_gosink = 0.1
    pmax = argmin(abs.(Re_gos .- 1)) # 63 is where Re_gosink = 1
    pmin = argmin(abs.(Re_gos .- 0.1)) # 15 is where Re_gosink = 0.01
    for indx in eachindex(vRHS)
        if Re_gos[indx] < Re_gos[pmin]
            vmixed[indx] = vLHS[indx]
        elseif Re_gos[indx] < Re_gos[pmax]
            fact = (log(Re_gos[indx]) - log(Re_gos[pmin]))^2 ./ ( (log(Re_gos[indx]) - log(Re_gos[pmin]))^2 .+ (log(Re_gos[pmax]) - log(Re_gos[indx]))^2)
            vmixed[indx] = vLHS[indx] * (1-fact) + vRHS[indx] * fact
        else
            vmixed[indx] = vRHS[indx]
        end
    end
    return vmixed
end

function find_Vi(R, H)
    return π*R^2*H
end

end