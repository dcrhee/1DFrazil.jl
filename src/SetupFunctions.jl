module SetupFunctions

export find_mean_collision_radius, find_collision_efficiency, find_β_α_Vol_constants, find_steady_velocity_via_iteration, find_Vi

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
            αs[Vindx] = findlast(Vpostcoll .- Vs .>= 0)
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

function find_steady_velocity_via_iteration(Rs, ν, ρₒ, ρᵢ, grav, aspect_ratio)
    vginitial = range(start = 1e-8, stop = 1e-1, step = 1e-9)
    vgs = zeros(length(Rs))
    for Rindx in eachindex(Rs)
        rᵢ = Rs[Rindx]
        Reₚ = 2 * rᵢ * vginitial/ν
        logCD = 1.386 .- 0.892*log10.(Reₚ) .+ 0.111 * (log10.(Reₚ)).^2
        CD = exp.(logCD)
        #print(vginitial.^2)
        #print(4*(ρₒ - ρᵢ)/ρₒ*grav * rᵢ/aspect_ratio * 1 ./CD)
        #print(maximum(4*(ρₒ - ρᵢ)/ρₒ * grav * rᵢ/aspect_ratio * 1 ./CD))
        #print(abs.(vginitial.^2 .- 4*(ρₒ - ρᵢ)/ρₒ*grav * rᵢ/aspect_ratio * 1 ./CD))
        vgindx = argmin(abs.(vginitial.^2 .- 4*(ρₒ - ρᵢ)/ρₒ*grav * rᵢ/aspect_ratio * 1 ./CD))
        vgs[Rindx] = vginitial[vgindx]
    end
    return vgs
end

function find_Vi(R, H)
    return π*R^2*H
end

end