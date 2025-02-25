using Oceananigans
using CairoMakie
using Oceananigans.Units: minute, minutes, hours
using Oceananigans.BuoyancyModels: g_Earth

using Oceananigans.AbstractOperations: ∂z
using Printf
using Statistics
using Oceananigans.AbstractOperations
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat, αₛ # these constants can be called inside any function

# 1. calculate the collision rate
# 2. calculate how the size distribution changes

# setup grid: choose 128 data points
depth = 0.20
grid = RectilinearGrid(size=128, z=(-depth, 0), topology=(Flat, Flat, Bounded))

# constants/parameters

# choose radius intervals
Rmin = 1e-5
Rmax = Rmin*10
Rstep = Rmin/10
aspect_ratio = 5

Rs = range(Rmin, stop=Rmax, step=Rstep)
Hs = 2*Rs/aspect_ratio
Vs = π * Rs.^2 .* Hs
Vs = range(1, stop=10, step=1)

V₁ = Vs[1]
Vₙ = Vs[end]

ϵ = 10^-3 # m²s⁻³
ν = 10^-6
nmax = 10^7
ζ = 1 # number of new crystals formed per collision
# work out the indices of the class to count crystal collisions from
Vrem = ζ * V₁

function find_β_α_Vol_constants()
    # matrix of which crystal size the crystal breaks off into
    αs = zeros(length(Vs))
    for Vindx in eachindex(Vs)
        Vpostcoll = Vs[Vindx] - Vrem
        print(Vpostcoll)
        if Vpostcoll >= V₁
            αs[Vindx] = findlast(Vpostcoll .- Vs .>= 0)
        end
    end

    # matrix of which crystal breaking introduces crystals of this size
    βs = zeros(length(Vs))
    for Vindx in eachindex(Vs)
        try
            βs[Vindx] = findfirst(x -> x == Vindx, αs)
        catch
        end
    end
    # coefficient for crystal of size class i breaking
    αVolconst = zeros(length(Vs))
    for Vindx in eachindex(αVolconst)
        try
            αVolconst[Vindx] = (1 + Vs[Vindx - αs[Vindx]]/(Vs[Vindx] - Vs[Vindx - αs[Vindx]]))
        catch
        end
    end
    αVolconst = zeros(length(Vs))
    for Vindx in eachindex(αVolconst)
        try
            αVolconst[Vindx] = -ζ * V₁/Vs[Vindx] * (1 + Vs[Vindx - αs[Vindx]]/(Vs[Vindx] - Vs[Vindx - αs[Vindx]]))
        catch
        end
    end
    # coefficient for crystal breaking into size class i
    βVolconst = zeros(length(Vs))
    for Vindx in eachindex(βVolconst)
        try
            βVolconst[Vindx] = ζ * V₁/(Vs[Vindx + βs[Vindx]] - Vs[Vindx])
        catch
        end
    end

    return αs, βs, αVolconst, βVolconst
end

αs, βs, αVolconst, βVolconst  = find_β_α_Vol_constants()

end_name = "no_vel"

collision_velocity_parameterisation_num = 1 # 1 is old cylinder, 2 is new cylinder, 3 is new spherical
concentration_parameterisation_num = 1 # 1 is mean n, 2 is sum over nj
crystal_size_collision_redistribution = 1 # 1 is old redistribution, 2 is new redistribution

#n1_boundary_conditions = FieldBoundaryConditions(bottom = FluxBoundaryCondition(nothing), top = FluxBoundaryCondition(nothing))
#n2_boundary_conditions = FieldBoundaryConditions(bottom = FluxBoundaryCondition(nothing), top = FluxBoundaryCondition(nothing))
#n3_boundary_conditions = FieldBoundaryConditions(bottom = FluxBoundaryCondition(nothing), top = FluxBoundaryCondition(nothing))

n1_boundary_conditions = FieldBoundaryConditions(top = GradientBoundaryCondition(0))
n2_boundary_conditions = FieldBoundaryConditions(top = GradientBoundaryCondition(0))
n3_boundary_conditions = FieldBoundaryConditions(top = GradientBoundaryCondition(0))

T_boundary_conditions = FieldBoundaryConditions(bottom = FluxBoundaryCondition(nothing), top = FluxBoundaryCondition(nothing))
S_boundary_conditions = FieldBoundaryConditions(bottom = FluxBoundaryCondition(nothing), top = FluxBoundaryCondition(nothing))

coriolis = FPlane(f=-1.4e-4) # s⁻¹


function find_Vi(R, H)
    return π*R^2*H
end

function find_density(T, S, ρ₀ = 1027.0, β = 0.0078, α = 1.67*10^(-4))
    ρ = ρ₀ * (1.0256550500000001 - α*T + β*S)
    return ρ
end

function find_steady_velocity(Rᵢ)
    uᵢ = 30*Rᵢ^(1.2)
    return uᵢ
end

function find_growth_rate(T, Rᵢ, H)
    G = kl*Nu/(ρᵢ*Lat) * (Tf - T) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
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
    Rᵢ = Rs[indx]
    H = Hs[indx]
    ρ = find_density(T, S)
    Tconstᵢ = kl*Nu/(ρ*cᴾ) * (Tf - T) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
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
        vcoll = sqrt(2/π)*(1/2*exp(-vᵣ^2/(2*vₜ^2)) + sqrt(π)/(2*sqrt(2))*(vᵣ/vₜ + vₜ/vᵣ)*erf(vᵣ/(sqrt(2)*vₜ)))
    end
    return vcoll
end

function nintermediate_forcing_func(i, j, k, grid, clock, model_fields, indx)
    ρ = find_density(model_fields.T, model_fields.S)
    uᵢ = find_steady_velocity(Rs[indx])
    ζᵢ = find_zeta(Rs[indx])

    # Compute derivative for nᵢ
    nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`

    # Extend velocity array (avoid index errors)
    uᵢ = find_steady_velocity(Rs[indx])
    udn_dz = -uᵢ /depth #use the faces below

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
    if concentration_parameterisation_num == 1
        Fenc = zeros(1, 1, 128)
        for Rindx in eachindex(Rs)
            nRindx = getfield(model_fields, Symbol("n$Rindx"))
            uRindx = find_steady_velocity(Rs[Rindx]) # find rise velocity of crystal j
            if collision_velocity_parameterisation_num == 2
                vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(Rs[indx] + Rs[Rindx]) # cylindrical approximation
            else
                vₜ = sqrt(ϵ/(15ν))*(Rs[indx] + Rs[Rindx]) # spherical approximation
            end
            vcoll = find_vcoll(uRindx - uᵢ, vₜ)
            if collision_velocity_parameterisation_num == 3
                Fenc .+= 2*π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
            else
                Fenc .+= π*(Rs[indx] + Rs[Rindx])^2*nRindx * vcoll
            end
        end    
    else
        Fenc = zeros(1, 1, 128)
        if collision_velocity_parameterisation_num == 2
            vₜ = sqrt(ϵ/(15ν)*(sqrt(π/2) + sqrt(2/π)))*(2*Rs[indx]) # cylindrical approximation
        else
            vₜ = sqrt(ϵ/(15ν))*(2*Rs[Rindx]) # spherical approximation
        end
        vcoll = find_vcoll(uᵢ, vₜ)
        nTotal = 0
        for Rindx in eachindex(Rs)
            nRindx = getfield(model_fields, Symbol("n$Rindx"))
            nTotal += sum(nRindx)
        end
        ntot = min(nTotal, nmax)

        if collision_velocity_parameterisation_num == 3
            Fenc = 2*π*(Rs[indx])^2 * vcoll * ntot
        else
            Fenc = π*(Rs[indx])^2 * vcoll * ntot
        end
    end
    Fcoll = Fend .* nᵢ

    # crystal size redistribution
    if crystal_size_collision_redistribution == 1
        dn_coll = - V₁/Vᵢ * Fcoll
    elseif Vᵢ > Vₙ - ζ*Vᵢ
        dn_coll = -ζ * V₁/Vᵢ * (1 + Vs[indx - αs[indx]]/(Vᵢ - Vs[indx - αs[indx]])) * Fcoll
    else
        dn_coll = -ζ * V₁/Vᵢ * (1 + Vs[indx - αs[indx]]/(Vᵢ - Vs[indx - αs[indx]])) * Fcoll + ζ * V₁/(Vs[indx + β] - Vᵢ) * (1 + Vs[indx - α]/(Vᵢ - Vs[indx - α])) * Fcollᵦ
    end

    # Apply logic for growth and melting
    if G_im1[i, j, k] > 0  # Growth case
        return @inbounds udn_dz - (G_ip1[i, j, k] * n_ip1[i, j, k] / (V_ip1 - Vᵢ) - G_im1[i, j, k] * n_im1[i, j, k] / (Vᵢ - V_im1)) + dn_coll #udn_dz[i, j, k] - (G_ip1[i, j, k] * n_ip1[i, j, k] / (V_ip1 - Vᵢ) - G_im1[i, j, k] * n_im1[i, j, k] / (Vᵢ - V_im1)) - ζᵢ * nᵢ * Fcoll
    else  # Melt case
        return @inbounds udn_dz - (Gᵢ[i, j, k] * nᵢ[i, j, k] / (V_ip1 - Vᵢ) - G_im1[i, j, k] * n_im1[i, j, k] / (Vᵢ - V_im1)) + dn_coll #udn_dz[i, j, k] - (Gᵢ[i, j, k] * nᵢ[i, j, k] / (V_ip1 - Vᵢ) - G_im1[i, j, k] * n_im1[i, j, k] / (Vᵢ - V_im1)) - ζᵢ * nᵢ * Fcoll
    end
end

function n1_forcing_func(i, j, k, grid, clock, model_fields, indx)
    uᵢ = find_steady_velocity(Rs[indx])
    udn_dz = -uᵢ /depth #use the faces below

    # growth term
    G₁ = find_growth_rate(model_fields.T, Rs[indx], Hs[indx])
    G₂ = find_growth_rate(model_fields.T, Rs[indx+1], Hs[indx+1])
    V₁ = find_Vi(Rs[indx], Hs[indx])
    V₂ = find_Vi(Rs[indx+1], Hs[indx+1])
    if G₁[i, j, k] > 0 # growth
        return @inbounds udn_dz - G₁[i, j, k]*model_fields.n1[i, j, k]/(V₂ - V₁) #@inbounds udn_dz[i, j, k] - G₁[i, j, k]*model_fields.n1[i, j, k]/(V₂ - V₁)
    else # melt
        return @inbounds udn_dz - (G₂[i, j, k]*model_fields.n2[i, j, k]/(V₂ - V₁) - G₁[i, j, k]*model_fields.n1[i, j, k]/V₁) #@inbounds udn_dz[i, j, k] - (G₂[i, j, k]*model_fields.n2[i, j, k]/(V₂ - V₁) - G₁[i, j, k]*model_fields.n1[i, j, k]/V₁)
    end
end
function n3_forcing_func(i, j, k, grid, clock, model_fields, indx)

    ρ = find_density(model_fields.T, model_fields.S)
    uᵢ = find_steady_velocity(Rs[indx])
    udn_dz = -uᵢ /depth #use the faces below

    # growth term
    G₂ = find_growth_rate(model_fields.T, Rs[indx-1], Hs[indx-1])
    G₃ = find_growth_rate(model_fields.T, Rs[indx], Hs[indx])
    V₂ = find_Vi(Rs[indx-1], Hs[indx-1])
    V₃ = find_Vi(Rs[indx], Hs[indx])
    if G₂[i, j, k] > 0 # growth
        return @inbounds udn_dz + G₂[i, j, k]*model_fields.n2[i, j, k]/(V₃ - V₂) #udn_dz[i, j, k] + G₂[i, j, k]*model_fields.n2[i, j, k]/(V₃ - V₂)
    else # melt
        return @inbounds udn_dz + (G₃[i, j, k]*model_fields.n3[i, j, k])/(V₃ - V₂) # udn_dz[i, j, k] + (G₃[i, j, k]*model_fields.n3[i, j, k])/(V₃ - V₂)
    end
end

function T_forcing_func(z, t, T, S)
    Tadd = zeros(1, 1, 128)
    for Rindx in eachindex(Rs)
        nRindx = getfield(model_fields, Symbol("n$Rindx"))
        Tconstᵢ = temperature_forcing_constant(T, S, indx)
        Tadd .+= nRindx*Tconstᵢ
    end
    return Tadd
end

function S_forcing_func(z, t, T, S, n1, n2, n3, p)
    Sconst₁ = salinity_forcing_constant(T, S, p.R₁)
    Sconst₂ = salinity_forcing_constant(T, S, p.R₂)
    Sconst₃ = salinity_forcing_constant(T, S, p.R₃)
    return Sconst₁ * n1 + + Sconst₂ * n2 + Sconst₃ * n3
end

n1_forcing = Forcing(n1_forcing_func, discrete_form=true, parameters = 1)
n2_forcing = Forcing(nintermediate_forcing_func, discrete_form=true, parameters = 2)
n3_forcing = Forcing(n3_forcing_func, discrete_form=true, parameters = 3)
T_forcing = Forcing(T_forcing_func, parameters=(cᴾ = cᴾ, k = kl, Nu = Nu, ρ=ρₒ, Tf = Tf), field_dependencies=(:T, :S, :n1, :n2, :n3))
S_forcing = Forcing(S_forcing_func, parameters=(cᴾ = cᴾ, α = α, k = kl, Nu = Nu, ρ=ρₒ, ρᵢ=ρᵢ, R₁ = Rs[1], R₂ = Rs[2], R₃ = Rs[3], H = 0.0004, Tf = Tf), field_dependencies=(:T, :S, :n1, :n2, :n3))

model = NonhydrostaticModel(; grid, coriolis,
advection = WENO(),
timestepper = :RungeKutta3,
tracers = (:T, :S, :n1, :n2, :n3),
buoyancy = SeawaterBuoyancy(),
closure = SmagorinskyLilly(Pr = 1, Cb = 1 / 1),
forcing=(n1=n1_forcing, n2=n2_forcing, n3=n3_forcing, S=S_forcing, T=T_forcing),
#stokes_drift = UniformStokesDrift(∂z_uˢ=∂z_uˢ),
boundary_conditions = (T=T_boundary_conditions, S=S_boundary_conditions, n1 = n1_boundary_conditions))


u, v, w = model.velocities
width = depth/10
nᵢ(z) = 1e6*exp(-(z+depth/2)^2 / (2width^2))

# set the initial conditions
Tᵢ(z) = Tf - 0.01 #* Ξₜ(z)
set!(model, w=0, T=Tᵢ, n1 = nᵢ, n2 = nᵢ, n3 = nᵢ, S=35)

simulation = Simulation(model, Δt=1.0, stop_time=2hours)

conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=1minute)


function progress(simulation)
    u, v, w = simulation.model.velocities

    # Print a progress message
    msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, wall time: %s\n",
    iteration(simulation),
    prettytime(time(simulation)),
    prettytime(simulation.Δt),
    maximum(abs, u), maximum(abs, v), maximum(abs, w),
    prettytime(simulation.run_wall_time))

    @info msg

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))

output_interval = 0.1minutes

fields_to_output = merge(model.velocities, model.tracers, (; νₑ=model.diffusivity_fields.νₑ))

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

run!(simulation)