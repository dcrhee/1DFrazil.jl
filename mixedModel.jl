using Oceananigans
using CairoMakie
using Oceananigans.Units: minute, minutes, hours
using SpecialFunctions
#using Oceananigans.BuoyancyModels: g_Earth

using Oceananigans.AbstractOperations: ∂z
using Printf
using Statistics
using Oceananigans.AbstractOperations
include("Constants.jl")
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat, αₛ, grav # these constants can be called inside any function

# setup grid: choose 128 data points
depth = 0.20
numz = 1 #128
#grid = RegularRectilinearGrid(size=(1, 1, 1), extent=(1.0, 1.0, 1.0))
grid = RectilinearGrid(size=1, z = (-1, 0), topology=(Flat, Flat, Periodic))
volume = 1

# run with the higher values of epsilon and see if it makes a difference using the different formulations

# constants/parameters

# choose radius intervals
aspect_ratio = 50

Rs = [0.01, 0.05, 0.15, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 2] .* 1e-3
Hs = 2*Rs/aspect_ratio
Vs = π * Rs.^2 .* Hs
#Vs = range(1, stop=10, step=1)

V₁ = Vs[1]
Vₙ = Vs[end]

ϵ = 1e-5#7.4 * 1e-6 #10^-3 # m²s⁻³
ν = 1.95 * 1e-6
Cᵢₙ =  4 * 1e-8
nmax = 10^20#10^3 * volume
ζ = 1 # number of new crystals formed per collision
# work out the indices of the class to count crystal collisions from
Vrem = ζ * V₁

function progress(simulation)
    u, v, w = simulation.model.velocities
    T = simulation.model.tracers.T
    n1 = simulation.model.tracers.n1
    n2 = simulation.model.tracers.n2
    n3 = simulation.model.tracers.n3
    n4 = simulation.model.tracers.n4
    n5 = simulation.model.tracers.n5
    n6 = simulation.model.tracers.n6
    n7 = simulation.model.tracers.n7
    n8 = simulation.model.tracers.n8
    n9 = simulation.model.tracers.n9
    n10 = simulation.model.tracers.n10

    # Print a progress message
    #msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, wall time: %s\n",
    msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, n1 = %.1e, n2 = %.1e, n3 = %.1e, n4 = %.1e, n5 = %.1e, n6 = %.1e, n7 = %.1e, n8 = %.1e, n9 = %.1e, n10 = %.1e, Tmin = %.5f, Tmax = %.5f, wall time: %s\n",
    iteration(simulation),
    prettytime(time(simulation)),
    prettytime(simulation.Δt),
    maximum(abs, u), maximum(abs, v), maximum(abs, w),
    minimum(n1), minimum(n2), minimum(n3), minimum(n4), minimum(n5), minimum(n6), minimum(n7), minimum(n8), minimum(n9), minimum(n10), 
    #maximum(n1), maximum(n2), maximum(n3), maximum(n4), maximum(n5), maximum(n6), maximum(n7), maximum(n8), maximum(n9), maximum(n10), 
    minimum(T), maximum(T),
    prettytime(simulation.run_wall_time))

    @info msg

    return nothing
end

function find_β_α_Vol_constants()
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

αs, βs, αVolconst, βVolconst  = find_β_α_Vol_constants()

function find_steady_velocity_via_iteration()
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

vgs = find_steady_velocity_via_iteration()

#n1_boundary_conditions = FieldBoundaryConditions(top = GradientBoundaryCondition(0))
#n2_boundary_conditions = FieldBoundaryConditions(top = GradientBoundaryCondition(0))
#n3_boundary_conditions = FieldBoundaryConditions(top = GradientBoundaryCondition(0))

#T_boundary_conditions = FieldBoundaryConditions(bottom = FluxBoundaryCondition(nothing), top = FluxBoundaryCondition(nothing))
#S_boundary_conditions = FieldBoundaryConditions(bottom = FluxBoundaryCondition(nothing), top = FluxBoundaryCondition(nothing))

coriolis = FPlane(f=-1.4e-4) # s⁻¹

# 1, 1
# 1, 2
# 2, 1
# 2, 2

cvelnum = 3
pnum = 2

#for cvelnum = [2, 3] #[1, 2, 3]
    #for pnum = [1, 2]
        end_name = "match_Feltham_no_n_max_eps_0.00001_depth_0.2"

        collision_velocity_parameterisation_num = cvelnum # 1 is old cylinder, 2 is new cylinder, 3 is new spherical
        concentration_parameterisation_num = pnum # 1 is mean n, 2 is sum over nj
        crystal_size_collision_redistribution = 1 # 1 is old redistribution, 2 is new redistribution
        effective_radius = true # add in their effective radius

        if collision_velocity_parameterisation_num == 2
            end_name = end_name * "_new_cyl"
        elseif collision_velocity_parameterisation_num == 3
            end_name = end_name * "_new_spherical"
        end
        if concentration_parameterisation_num == 2
            end_name = end_name * "_sum_nj"
        end

        function find_Vi(R, H)
            return π*R^2*H
        end
        
        function find_density(T, S)
            ρ = ρₒ * (1 + 7.86 * 1e-4 * (S - 34.5) - 3.87 * 1e-5 * (T + 2) )
            return ρ
        end
        
        function find_steady_velocity(indx)
            #uᵢ = 30*Rᵢ^(1.2)
            uᵢ = vgs[indx]
            return uᵢ
        end
        
        function find_growth_rate(T, Rᵢ, H)
            #G = kl*Nu/(ρᵢ*Lat) * (Tf - T) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
            G = kl*Nu/(ρᵢ*Lat) * (Tf - T) * 2π * H
            return G
        end
        
        function find_zeta(Rindx)
            # geometric factor to account for the different volumes of each size class as you move from one to the other
            V₁ = find_Vi(Rs[1], Hs[1])
            Vᵢ = find_Vi(Rs[Rindx], Hs[Rindx])
            ζᵢ = V₁/Vᵢ
            return ζᵢ
        end
        
        function temperature_forcing_constant(T, S, indx)
            # constant in front of each concentration
            #Rᵢ = Rs[indx]
            H = Hs[indx]
            ρ = find_density(T, S)
            Tconstᵢ = kl*Nu/(volume * ρ*cᴾ) * (Tf - T) * 2π * H
            #Tconstᵢ = kl*Nu/(volume * ρ*cᴾ) * (Tf - T) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
            return Tconstᵢ
        end
        
        function salinity_forcing_constant(T, S, indx)
            # constant in front of each concentration
            Rᵢ = Rs[indx]
            H = Hs[indx]
            ρ = find_density(T, S)
            Sconstᵢ = S*(1-αₛ)*kl*Nu/(ρ*Lat) * (Tf - T) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
            return Sconstᵢ
        end
        
        function find_vcoll(vᵣ, vₜ)
            if collision_velocity_parameterisation_num == 1
                vcoll = sqrt(vᵣ^2 + vₜ^2)
            elseif collision_velocity_parameterisation_num == 2
                vcoll = sqrt(vᵣ^2 + (π+2)^2/(2π)*vₜ^2)
            else
                if vᵣ == 0
                    vcoll = sqrt(2/π)*vₜ
                else
                    vcoll = sqrt(2/π)*vₜ*(1/2*exp(-vᵣ^2/(2*vₜ^2)) + sqrt(π)/(2*sqrt(2))*(vᵣ/vₜ + vₜ/vᵣ)*erf(vᵣ/(sqrt(2)*vₜ)))
                end
            end
            return vcoll
        end
        
        function nintermediate_forcing_func(i, j, k, grid, clock, model_fields, indx)
            ρ = find_density(model_fields.T, model_fields.S)
            uᵢ = find_steady_velocity(indx)
            #ζᵢ = find_zeta(Rs[indx])
        
            # Compute derivative for nᵢ
            nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
            nᵢ = max.(nᵢ, 0)
        
            # Extend velocity array (avoid index errors)
            udn_dz = -nᵢ .* uᵢ /depth #use the faces below
        
            # Access nᵢ₋₁ and nᵢ₊₁ safely
            n_im1 = getfield(model_fields, Symbol("n$(indx-1)"))  # Field nᵢ₋₁
            n_ip1 = getfield(model_fields, Symbol("n$(indx+1)"))  # Field nᵢ₊₁
        
            # Compute growth rates and Vi values dynamically
            Gᵢ = find_growth_rate(model_fields.T, Rs[indx], Hs[indx])
            G_im1 = find_growth_rate(model_fields.T, Rs[indx-1], Hs[indx-1])
            G_ip1 = find_growth_rate(model_fields.T, Rs[indx+1], Hs[indx+1])
        
            Vᵢ = find_Vi(Rs[indx], Hs[indx])
            V_im1 = find_Vi(Rs[indx-1], Hs[indx-1])
            V_ip1 = find_Vi(Rs[indx+1], Hs[indx+1])
        
            # get collision frquency
            if αs[indx] == 0
                Fcoll = zeros(1, 1, numz)
            else
                if concentration_parameterisation_num == 2
                    Fenc = zeros(1, 1, numz)
                    for Rindx in eachindex(Rs)
                        nRindx = getfield(model_fields, Symbol("n$Rindx"))
                        uRindx = find_steady_velocity(Rindx) # find rise velocity of crystal j
                        if collision_velocity_parameterisation_num == 2
                            if effective_radius
                                vₜ = (3/(2*aspect_ratio))^(1/3)*sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[indx] + Rs[Rindx]) # cylindrical approximation
                            else
                                vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[indx] + Rs[Rindx]) # cylindrical approximation
                            end
                        else
                            if effective_radius
                                vₜ = (3/(2*aspect_ratio))^(1/3)*sqrt(ϵ/(15ν))*(Rs[indx] + Rs[Rindx]) # spherical approximation
                            else
                                vₜ = sqrt(ϵ/(15ν))*(Rs[indx] + Rs[Rindx]) # spherical approximation
                            end
                        end
                        vcoll = find_vcoll(uRindx - uᵢ, vₜ)
                        if collision_velocity_parameterisation_num == 3
                            if effective_radius
                                Fenc .+= (3/(2*aspect_ratio))^(2/3)*2*π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
                            else
                                Fenc .+= 2*π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
                            end
                        else
                            if effective_radius 
                                Fenc .+= (3/(2*aspect_ratio))^(2/3)*π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
                            else
                                Fenc .+= π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
                            end
                        end
                    end    
                else
                    Fenc = zeros(1, 1, numz)
                    if collision_velocity_parameterisation_num == 2
                        if effective_radius
                            vₜ = (3/(2*aspect_ratio))^(1/3)*sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[indx]) # cylindrical approximation
                        else
                            vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[indx]) # cylindrical approximation
                        end
                    else
                        if effective_radius
                            vₜ = (3/(2*aspect_ratio))^(1/3)*sqrt(ϵ/(15ν))*(2*Rs[indx]) # spherical approximation
                        else
                            vₜ = sqrt(ϵ/(15ν))*(2*Rs[indx]) # spherical approximation
                        end
                    end
                    vcoll = find_vcoll(uᵢ, vₜ)
                    nTotal = zeros(1, 1, numz)
                    for Rindx in eachindex(Rs)
                        nRindx = getfield(model_fields, Symbol("n$Rindx"))
                        nTotal .+= nRindx
                    end
                    ntot = min.(nTotal, nmax)
        
                    if collision_velocity_parameterisation_num == 3
                        if effective_radius
                            Fenc =  (3/(2*aspect_ratio))^(2/3)*2*π*(Rs[indx])^2 * vcoll * ntot
                        else
                            Fenc = 2*π*(Rs[indx])^2 * vcoll * ntot
                        end
                    else
                        if effective_radius
                            Fenc = (3/(2*aspect_ratio))^(2/3)*π*(Rs[indx])^2 * vcoll * ntot
                        else
                            Fenc = π*(Rs[indx])^2 * vcoll * ntot
                        end
                    end
                end
                Fcoll = Fenc .* nᵢ
            end
        
            # find Fcollᵦ
            if crystal_size_collision_redistribution == 2
                β = βs[indx]
                if β == 0
                    Fcollᵦ == zeros(1, 1, numz)
                else
                    uᵦ = find_steady_velocity(β)
                    
                    nᵦ = getfield(model_fields, Symbol("n$β"))  # Dynamically get field `nᵢ`
                
                    
                    if concentration_parameterisation_num == 2
                        Fenc = zeros(1, 1, numz)
                        for Rindx in eachindex(Rs)
                            nRindx = getfield(model_fields, Symbol("n$Rindx"))
                            uRindx = find_steady_velocity(Rindx) # find rise velocity of crystal j
                            if collision_velocity_parameterisation_num == 2
                                if effective_radius
                                    vₜ = (3/(2*aspect_ratio))^(1/3)*sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[β] + Rs[Rindx]) # cylindrical approximation
                                else
                                    vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[β] + Rs[Rindx]) # cylindrical approximation
                                end
                            else
                                if effective_radius
                                    vₜ = (3/(2*aspect_ratio))^(1/3)*sqrt(ϵ/(15ν))*(Rs[β] + Rs[Rindx]) # spherical approximation
                                else
                                    vₜ = sqrt(ϵ/(15ν))*(Rs[β] + Rs[Rindx]) # spherical approximation
                                end
                            end
                            vcoll = find_vcoll(uRindx - uᵦ, vₜ)
                            if collision_velocity_parameterisation_num == 3
                                if effective_radius
                                    Fenc .+= (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[β] + Rs[Rindx])^2*nRindx * vcoll
                                else
                                    Fenc .+= 2*π*(Rs[β] + Rs[Rindx])^2*nRindx * vcoll
                                end
                            else
                                if effective_radius
                                    Fenc .+= (3/(2*aspect_ratio))^(2/3) * π*(Rs[β] + Rs[Rindx])^2*nRindx * vcoll
                                else
                                    Fenc .+= π*(Rs[β] + Rs[Rindx])^2*nRindx * vcoll
                                end
                            end
                        end    
                    else
                        Fenc = zeros(1, 1, numz)
                        if collision_velocity_parameterisation_num == 2
                            if effective_radius
                                vₜ = (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[β]) # cylindrical approximation
                            else
                                vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[β]) # cylindrical approximation
                            end
                        else
                            if effective_radius
                                vₜ = (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν))*(2*Rs[β]) # spherical approximation
                            else
                                vₜ = sqrt(ϵ/(15ν))*(2*Rs[β]) # spherical approximation
                            end
                        end
                        vcoll = find_vcoll(uᵦ, vₜ)
                        nTotal = zeros(1, 1, numz)
                        for Rindx in eachindex(Rs)
                            nRindx = getfield(model_fields, Symbol("n$Rindx"))
                            nTotal .+= nRindx
                        end
                        ntot = min.(nTotal, nmax)
        
                        if collision_velocity_parameterisation_num == 3
                            if effective_radius
                                Fencᵦ = (3/(2*aspect_ratio))^(2/3) *  2*π*(Rs[β])^2 * vcoll * ntot
                            else
                                Fencᵦ = 2*π*(Rs[β])^2 * vcoll * ntot
                            end
                        else
                            if effective_radius
                                Fencᵦ = (3/(2*aspect_ratio))^(2/3) *  π*(Rs[β])^2 * vcoll * ntot
                            else
                                Fencᵦ = π*(Rs[β])^2 * vcoll * ntot
                            end
                        end
                    end
                    Fcollᵦ = Fencᵦ .* nᵦ
                end
            end
            
        
        
            # crystal size redistribution
            if crystal_size_collision_redistribution == 1
                dn_coll = - V₁/Vᵢ * Fcoll/volume
            else
                dn_coll = (αVolconst[indx] * Fcoll .+ βVolconst[indx] * Fcollᵦ)/volume
            end
            #print("n2", maximum(dn_coll))
        
            # Apply logic for growth and melting
            if G_im1[i, j, k] > 0  # Growth case
                final_conc = - (Gᵢ[i, j, k] * nᵢ[i, j, k] / (V_ip1 - Vᵢ) - G_im1[i, j, k] * n_im1[i, j, k] / (Vᵢ - V_im1)) + dn_coll[i, j, k] + udn_dz[i, j, k]
                return @inbounds final_conc
            else  # Melt case
                final_conc = - (G_ip1[i, j, k] * n_ip1[i, j, k] / (V_ip1 - Vᵢ) - Gᵢ[i, j, k] * nᵢ[i, j, k] / (Vᵢ - V_im1)) + dn_coll[i, j, k] + udn_dz[i, j, k]
                return @inbounds final_conc
            end
        end
        
        function n1_forcing_func(i, j, k, grid, clock, model_fields, indx)
            uᵢ = find_steady_velocity(indx)
            nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
            nᵢ = max.(nᵢ, 0)
            udn_dz = -nᵢ .* uᵢ /depth #use the faces below
        
            # growth term
            G₁ = find_growth_rate(model_fields.T, Rs[indx], Hs[indx])
            G₂ = find_growth_rate(model_fields.T, Rs[indx+1], Hs[indx+1])
            V₁ = find_Vi(Rs[indx], Hs[indx])
            V₂ = find_Vi(Rs[indx+1], Hs[indx+1])
        
            # get collision frquency
            if concentration_parameterisation_num == 2
                Fenc = zeros(1, 1, numz)
                Fcollsum =  zeros(1, 1, numz)
                for jindx in 2:lastindex(Rs)
                    uⱼ = find_steady_velocity(jindx)
                    nⱼ = getfield(model_fields, Symbol("n$jindx"))  # Dynamically get field `nᵢ`
                    for Rindx in eachindex(Rs)
                        nRindx = getfield(model_fields, Symbol("n$Rindx"))
                        uRindx = find_steady_velocity(Rindx) # find rise velocity of crystal j
                        if collision_velocity_parameterisation_num == 2
                            if effective_radius 
                                vₜ = (3/(2*aspect_ratio))^(1/3) *  sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[jindx] + Rs[Rindx]) # cylindrical approximation
                            else
                                vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[jindx] + Rs[Rindx]) # cylindrical approximation
                            end
                        else
                            if effective_radius
                                vₜ = (3/(2*aspect_ratio))^(1/3) *  sqrt(ϵ/(15ν))*(Rs[jindx] + Rs[Rindx]) # spherical approximation
                            else
                                vₜ = sqrt(ϵ/(15ν))*(Rs[jindx] + Rs[Rindx]) # spherical approximation
                            end
                        end
                        vcoll = find_vcoll(uRindx - uⱼ, vₜ)
                        if collision_velocity_parameterisation_num == 3
                            if effective_radius
                                Fenc .+=  (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
                            else
                                Fenc .+= 2*π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
                            end
                        else
                            if effective_radius
                                Fenc .+= (3/(2*aspect_ratio))^(2/3) *  π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
                            else
                                Fenc .+= π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
                            end
                        end
                    end  
                    Fcollsum = Fcollsum .+ Fenc .* nⱼ
                end  
            else
                Fenc = zeros(1, 1, numz)
                Fcollsum =  zeros(1, 1, numz)
                for jindx in 2:lastindex(Rs)
                    uⱼ = find_steady_velocity(jindx)
                    nⱼ = getfield(model_fields, Symbol("n$jindx"))  # Dynamically get field `nᵢ`
                    if collision_velocity_parameterisation_num == 2
                        if effective_radius
                            vₜ = (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[jindx]) # cylindrical approximation
                        else
                            vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[jindx]) # cylindrical approximation
                        end
                    else
                        if effective_radius
                            vₜ = (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν))*(2*Rs[jindx]) # spherical approximation
                        else
                            vₜ = sqrt(ϵ/(15ν))*(2*Rs[jindx]) # spherical approximation
                        end
                    end
                    vcoll = find_vcoll(uⱼ, vₜ)
                    nTotal = zeros(1, 1, numz)
                    for Rindx in eachindex(Rs)
                        nRindx = getfield(model_fields, Symbol("n$Rindx"))
                        nTotal .+= nRindx
                    end
                    ntot = min.(nTotal, nmax)
        
                    if collision_velocity_parameterisation_num == 3
                        if effective_radius
                            Fenc = (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[jindx])^2 * vcoll * ntot
                        else
                            Fenc = 2*π*(Rs[jindx])^2 * vcoll * ntot
                        end
                    else
                        if effective_radius
                            Fenc =  (3/(2*aspect_ratio))^(2/3) *  π*(Rs[jindx])^2 * vcoll * ntot
                        else
                            Fenc = π*(Rs[jindx])^2 * vcoll * ntot
                        end
                    end
                    Fcollsum = Fcollsum .+ Fenc .* nⱼ
                end
            end
        
            
            if crystal_size_collision_redistribution == 1
                dn_coll = Fcollsum / volume
            else
                dn_coll = ζ * Fcollsum / volume
            end
            #print("n1", maximum(dn_coll))
        
            if G₁[i, j, k] > 0 # growth
                final_conc =   - G₁[i, j, k]*model_fields.n1[i, j, k]/(V₂ - V₁) + dn_coll[i, j, k] + udn_dz[i, j, k]
                return @inbounds final_conc
            else # melt
                final_conc = - (G₂[i, j, k]*model_fields.n2[i, j, k]/(V₂ - V₁) - G₁[i, j, k]*model_fields.n1[i, j, k]/V₁) + dn_coll[i, j, k] + udn_dz[i, j, k]
                
                return @inbounds  final_conc
            end
        end
        
        
        function nend_forcing_func(i, j, k, grid, clock, model_fields, indx)
        
            ρ = find_density(model_fields.T, model_fields.S)
            uᵢ = find_steady_velocity(indx)
            nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
            nᵢ = max.(nᵢ, 0)

            udn_dz = -nᵢ .* uᵢ /depth #use the faces below
        
            # growth term
            G₂ = find_growth_rate(model_fields.T, Rs[indx-1], Hs[indx-1])
            G₃ = find_growth_rate(model_fields.T, Rs[indx], Hs[indx])
            V₂ = find_Vi(Rs[indx-1], Hs[indx-1])
            Vᵢ = find_Vi(Rs[indx], Hs[indx])
        
            # get collision frquency
            if concentration_parameterisation_num == 2
                Fenc = zeros(1, 1, numz)
                for Rindx in eachindex(Rs)
                    nRindx = getfield(model_fields, Symbol("n$Rindx"))
                    uRindx = find_steady_velocity(Rindx) # find rise velocity of crystal j
                    if collision_velocity_parameterisation_num == 2
                        if effective_radius
                            vₜ = (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[indx] + Rs[Rindx]) # cylindrical approximation
                        else
                            vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[indx] + Rs[Rindx]) # cylindrical approximation
                        end
                    else
                        if effective_radius
                            vₜ =  (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν))*(Rs[indx] + Rs[Rindx]) # spherical approximation
                        else
                            vₜ = sqrt(ϵ/(15ν))*(Rs[indx] + Rs[Rindx]) # spherical approximation
                        end
                    end
                    vcoll = find_vcoll(uRindx - uᵢ, vₜ)
                    if collision_velocity_parameterisation_num == 3
                        if effective_radius
                            Fenc .+= (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll
                        else
                            Fenc .+= 2*π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll
                        end
                    else
                        if effective_radius
                            Fenc .+= (3/(2*aspect_ratio))^(2/3) * π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll
                        else
                            Fenc .+= π*(Rs[indx] + Rs[Rindx])^2 * nRindx * vcoll
                        end
                    end
                end    
            else
                Fenc = zeros(1, 1, numz)
                if collision_velocity_parameterisation_num == 2
                    if effective_radius
                        vₜ = (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[indx]) # cylindrical approximation
                    else
                        vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[indx]) # cylindrical approximation
                    end
                else
                    if effective_radius
                        vₜ = (3/(2*aspect_ratio))^(1/3) * sqrt(ϵ/(15ν))*(2*Rs[indx]) # spherical approximation
                    else
                        vₜ = sqrt(ϵ/(15ν))*(2*Rs[indx]) # spherical approximation
                    end
                end
                vcoll = find_vcoll(uᵢ, vₜ)
                nTotal = zeros(1, 1, numz)
                for Rindx in eachindex(Rs)
                    nRindx = getfield(model_fields, Symbol("n$Rindx"))
                    nTotal .+= nRindx
                end
                ntot = min.(nTotal, nmax)
        
                if collision_velocity_parameterisation_num == 3
                    if effective_radius
                        Fenc =  (3/(2*aspect_ratio))^(2/3) * 2*π*(Rs[indx])^2 * vcoll * ntot
                    else
                        Fenc = 2*π*(Rs[indx])^2 * vcoll * ntot
                    end
                else
                    if effective_radius
                        Fenc = (3/(2*aspect_ratio))^(2/3) *π*(Rs[indx])^2 * vcoll * ntot
                    else
                        Fenc = π*(Rs[indx])^2 * vcoll * ntot
                    end
                end
            end
            nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
            Fcoll = Fenc .* nᵢ
        
            # crystal size redistribution
            if crystal_size_collision_redistribution == 1
                dn_coll = - V₁/Vᵢ * Fcoll /volume
            else
                dn_coll = αVolconst[indx] * Fcoll /volume
            end
        
            if G₂[i, j, k] > 0 # growth
                final_conc = G₂[i, j, k]*model_fields.n9[i, j, k]/(Vᵢ - V₂) + dn_coll[i, j, k] + udn_dz[i, j, k]
                #print(", udn/dz= ", udn_dz[1, 1, 1], ", ")
                return @inbounds final_conc
            else # melt
                final_conc = (G₃[i, j, k]*model_fields.n10[i, j, k])/(Vᵢ - V₂) + dn_coll[i, j, k] + udn_dz[i, j, k]
                return @inbounds  final_conc
            end
        
        end
        
        function T_forcing_func(z, t, T, S, n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈, n₉, n10)
            Tconst₁ = temperature_forcing_constant(T, S, 1)
            Tconst₂ = temperature_forcing_constant(T, S, 2)
            Tconst₃ = temperature_forcing_constant(T, S, 3)
            Tconst4 = temperature_forcing_constant(T, S, 4)
            Tconst5 = temperature_forcing_constant(T, S, 5)
            Tconst6 = temperature_forcing_constant(T, S, 6)
            Tconst7 = temperature_forcing_constant(T, S, 7)
            Tconst8 = temperature_forcing_constant(T, S, 8)
            Tconst9 = temperature_forcing_constant(T, S, 9)
            Tconst10 = temperature_forcing_constant(T, S, 10)
            return Tconst₁ * n₁ + Tconst₂ * n₂ + Tconst₃ * n₃ + Tconst4 * n₄ + Tconst5 * n₅ + Tconst6 * n₆ + Tconst7 * n₇ + Tconst8 * n₈ + Tconst9 * n₉ + Tconst10 * n10
        end
        
        function S_forcing_func(z, t, T, S, n1, n2, n3)
            Sconst₁ = salinity_forcing_constant(T, S, 1)
            Sconst₂ = salinity_forcing_constant(T, S, 2)
            Sconst₃ = salinity_forcing_constant(T, S, 3)
            return Sconst₁ * n1 + + Sconst₂ * n2 + Sconst₃ * n3
        end

        n1_forcing = Forcing(n1_forcing_func, discrete_form=true, parameters = 1)
        n2_forcing = Forcing(nintermediate_forcing_func, discrete_form=true, parameters = 2)
        n3_forcing = Forcing(nintermediate_forcing_func, discrete_form=true, parameters = 3)
        n4_forcing = Forcing(nintermediate_forcing_func, discrete_form=true, parameters = 4)
        n5_forcing = Forcing(nintermediate_forcing_func, discrete_form=true, parameters = 5)
        n6_forcing = Forcing(nintermediate_forcing_func, discrete_form=true, parameters = 6)
        n7_forcing = Forcing(nintermediate_forcing_func, discrete_form=true, parameters = 7)
        n8_forcing = Forcing(nintermediate_forcing_func, discrete_form=true, parameters = 8)
        n9_forcing = Forcing(nintermediate_forcing_func, discrete_form=true, parameters = 9)
        n10_forcing = Forcing(nend_forcing_func, discrete_form=true, parameters = 10)
        T_forcing = Forcing(T_forcing_func, field_dependencies=(:T, :S, :n1, :n2, :n3, :n4, :n5, :n6, :n7, :n8, :n9, :n10))
        #S_forcing = Forcing(S_forcing_func, field_dependencies=(:T, :S, :n1, :n2, :n3))

        model = NonhydrostaticModel(; grid,
        advection = Centered(), #WENO(),
        timestepper = :RungeKutta3,
        buoyancy = SeawaterBuoyancy(),
        #closure = SmagorinskyLilly(Pr = 1, Cb = 1 / 1),
        tracers = (:T, :S, :n1, :n2, :n3, :n4, :n5, :n6, :n7, :n8, :n9, :n10),
        #buoyancy = SeawaterBuoyancy(),
        forcing=(n1=n1_forcing, n2=n2_forcing, n3=n3_forcing, n4=n4_forcing, n5=n5_forcing, n6=n6_forcing, n7=n7_forcing, n8=n8_forcing, n9=n9_forcing, n10=n10_forcing, T=T_forcing)) #,
        #boundary_conditions = (T=T_boundary_conditions, S=S_boundary_conditions, n1 = n1_boundary_conditions))

        u, v, w = model.velocities
        nᵢₙ = Cᵢₙ * aspect_ratio * volume / (2 * π * length(Rs)) * 1 ./ Rs.^3

        # set the initial conditions
        Tᵢ = Tf - 1e-4 #* Ξₜ(z)
        set!(model, u=0, v=0, w=0, T=Tᵢ, n1 = nᵢₙ[1], n2 = nᵢₙ[2], n3 = nᵢₙ[3], n4 = nᵢₙ[4], n5 = nᵢₙ[5], n6 = nᵢₙ[6], n7 = nᵢₙ[7], n8 = nᵢₙ[8], n9 = nᵢₙ[9], n10 = nᵢₙ[10], S=34.5)

        simulation = Simulation(model, Δt=1.0, stop_time=72hours)

        # Define the enforce_nonnegative_tracer function
        function enforce_nonnegative_tracer(simulation)
            n1_data = simulation.model.tracers.n1.data
            @inbounds n1_data .= max.(n1_data, 0)

            n2_data = simulation.model.tracers.n2.data
            @inbounds n2_data .= max.(n2_data, 0)

            n3_data = simulation.model.tracers.n3.data
            @inbounds n3_data .= max.(n3_data, 0)

            n4_data = simulation.model.tracers.n4.data
            @inbounds n4_data .= max.(n4_data, 0)

            n5_data = simulation.model.tracers.n5.data
            @inbounds n5_data .= max.(n5_data, 0)

            n6_data = simulation.model.tracers.n6.data
            @inbounds n6_data .= max.(n6_data, 0)

            n7_data = simulation.model.tracers.n7.data
            @inbounds n7_data .= max.(n7_data, 0)

            n8_data = simulation.model.tracers.n8.data
            @inbounds n8_data .= max.(n8_data, 0)

            n9_data = simulation.model.tracers.n9.data
            @inbounds n9_data .= max.(n9_data, 0)

            n10_data = simulation.model.tracers.n10.data
            @inbounds n10_data .= max.(n10_data, 0)

            return nothing
        end

        # Define a callback to enforce non-negative tracer values
        #function enforce_nonnegative_tracer(simulation)
            #@inbounds simulation.model.tracers.n1 .= max!.(simulation.model.tracers.n1, 0)  
            #@inbounds set!(simulation.model.tracers.n1) do n1
            #return max(n1, 0)
        #    set!(model, u=u, v=v, w=w, T=T, n1 = max(0, n1), n2 = n2, n3 = n3, n4 = n4, n5 = n5, n6 = n6, n7 = n7, n8 = n8, n9 = n9, n10 = n10, S=34.5)
            #end 
        #end

        # Create the simulation
        grid = RectilinearGrid(size=(32, 32, 32), extent=(1, 1, 1))
        
        # Add the callback to enforce non-negativity
        add_callback!(simulation, enforce_nonnegative_tracer, IterationInterval(1))
        add_callback!(simulation, progress, IterationInterval(20))

        conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=1minute)#0.01minute)

        #simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))

        output_interval = 1minutes

        fields_to_output = merge(model.velocities, model.tracers)

        simulation.output_writers[:fields] =
            JLD2OutputWriter(model, fields_to_output,
                            schedule = TimeInterval(output_interval),
                            filename = "1D_fields" * end_name * ".jld2",
                            overwrite_existing = true)

        u, v, w = model.velocities
        T = model.tracers.T
        S = model.tracers.S
        n1 = model.tracers.n1
        n2 = model.tracers.n2
        n3 = model.tracers.n3
        n4 = model.tracers.n4
        n5 = model.tracers.n5
        n6 = model.tracers.n6
        n7 = model.tracers.n7
        n8 = model.tracers.n8
        n9 = model.tracers.n9
        n10 = model.tracers.n10

        run!(simulation)
    #end
#end