using Oceananigans
using CairoMakie
using Oceananigans.Units: minute, minutes, hours
using SpecialFunctions
using Oceananigans.Fields: interpolate!
#using Oceananigans.BuoyancyModels: g_Earth

using Oceananigans.AbstractOperations: ∂z
using Printf
using Statistics
using Oceananigans.AbstractOperations
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat, αₛ, grav # these constants can be called inside any function

# setup grid: choose 128 data points
depth = 10
numz = 128 #128
#grid = RegularRectilinearGrid(size=(1, 1, 1), extent=(1.0, 1.0, 1.0))
grid = RectilinearGrid(size=numz, z = (-depth, 0), topology=(Flat, Flat, Bounded))
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

n_bcs = FieldBoundaryConditions(bottom=FluxBoundaryCondition(0), 
                                top=FluxBoundaryCondition(0))

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
    msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, n1 = (%.1e, %.1e), n2 = (%.1e, %.1e), n3 = (%.1e, %.1e), n4 = (%.1e, %.1e), n5 = (%.1e, %.1e), n6 = (%.1e, %.1e), n7 = (%.1e, %.1e), n8 = (%.1e, %.1e), n9 = (%.1e, %.1e), n10 = (%.1e, %.1e), Tmin = %.5f, Tmax = %.5f, wall time: %s\n",
    iteration(simulation),
    prettytime(time(simulation)),
    prettytime(simulation.Δt),
    maximum(abs, u), maximum(abs, v), maximum(abs, w),
    minimum(n1), maximum(n1), minimum(n2), maximum(n2), minimum(n3), maximum(n3), minimum(n4), maximum(n4), minimum(n5), maximum(n5), minimum(n6), maximum(n6), minimum(n7), maximum(n7), minimum(n8), maximum(n8), minimum(n9), maximum(n9), minimum(n10), maximum(n10), 
    #, , maximum(n3), maximum(n4), maximum(n5), maximum(n6), maximum(n7), maximum(n8), maximum(n9), maximum(n10), 
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

cvelnum = 1
pnum = 1

end_name = "match_Feltham_just_rise_noflux"

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


function nintermediate_forcing_func(i, j, k, grid, clock, model_fields, indx)
    ρ = find_density(model_fields.T, model_fields.S)
    #uᵢ = find_steady_velocity(indx)
    #ζᵢ = find_zeta(Rs[indx])

    # Compute derivative for nᵢ
    nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
    dn_dz = ∂z(nᵢ)
    #nᵢ = max.(nᵢ, 0)

    uᵢ_faces = Field{Face, Face, Face}(grid)  # Create a field on cell faces
    uᵢ = Field{Center, Center, Center}(grid)
    set!(uᵢ, find_steady_velocity(indx))
    interpolate!(uᵢ_faces, uᵢ)  
    udn_dz = -uᵢ_faces .* dn_dz #use the faces below
    # Enforce no-flux: Set dn_dz to zero at the boundaries
    if k == 1 || k == size(udn_dz, 3)  || k == size(udn_dz, 3) - 1  # Bottom or top boundary
        udn_dz[i, j, k] = 0  # Force ∂z(n) = 0 at boundaries
    end
    final_conc = udn_dz[i, j, k]
    return @inbounds final_conc
end

function n1_forcing_func(i, j, k, grid, clock, model_fields, indx)
    #uᵢ = find_steady_velocity(indx)
    nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
    dn_dz = ∂z(nᵢ)

    uᵢ_faces = Field{Face, Face, Face}(grid)  # Create a field on cell faces
    uᵢ = Field{Center, Center, Center}(grid)
    set!(uᵢ, find_steady_velocity(indx))
    interpolate!(uᵢ_faces, uᵢ)  
    udn_dz = -uᵢ_faces .* dn_dz #use the faces below
    # Enforce no-flux: Set dn_dz to zero at the boundaries
    if k == 1 || k == size(udn_dz, 3)  || k == size(udn_dz, 3) - 1 # Bottom or top boundary
        udn_dz[i, j, k] = 0  # Force ∂z(n) = 0 at boundaries
    end
    #print("u_i", uᵢ)
    #print("udn_dz", udn_dz)
    return @inbounds udn_dz[i, j, k]
    
end


function nend_forcing_func(i, j, k, grid, clock, model_fields, indx)

    ρ = find_density(model_fields.T, model_fields.S)
    #uᵢ = find_steady_velocity(indx)
    nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
    dn_dz = ∂z(nᵢ)
    uᵢ_faces = Field{Face, Face, Face}(grid)  # Create a field on cell faces
    uᵢ = Field{Center, Center, Center}(grid)
    set!(uᵢ, find_steady_velocity(indx))
    interpolate!(uᵢ_faces, uᵢ)  
    udn_dz = -uᵢ_faces .* dn_dz #use the faces below
    # Enforce no-flux: Set dn_dz to zero at the boundaries
    if k == 1 || k == size(udn_dz, 3) || k == size(udn_dz, 3) - 1  # Bottom or top boundary
        udn_dz[i, j, k] = 0  # Force ∂z(n) = 0 at boundaries
    end
    return @inbounds  udn_dz[i, j, k]
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

model = NonhydrostaticModel(; grid,
advection = Centered(), #WENO(),
timestepper = :RungeKutta3,
buoyancy = SeawaterBuoyancy(),
#closure = SmagorinskyLilly(Pr = 1, Cb = 1 / 1),
tracers = (:T, :S, :n1, :n2, :n3, :n4, :n5, :n6, :n7, :n8, :n9, :n10),
#buoyancy = SeawaterBuoyancy(),
forcing=(n1=n1_forcing, n2=n2_forcing, n3=n3_forcing, n4=n4_forcing, n5=n5_forcing, n6=n6_forcing, n7=n7_forcing, n8=n8_forcing, n9=n9_forcing, n10=n10_forcing),
boundary_conditions = (n1 = n_bcs, n2 = n_bcs, n3 = n_bcs, n4 = n_bcs, n5 = n_bcs, n6 = n_bcs, n7 = n_bcs, n8 = n_bcs, n9 = n_bcs, n10 = n_bcs))

u, v, w = model.velocities
#nᵢₙ = Cᵢₙ * aspect_ratio * volume / (2 * π * length(Rs)) * 1 ./ Rs.^3

width = depth/10
n(z) = 1e10*exp(-(z+depth/2)^2 / (2width^2))

Ξₜ(z) = randn()  # noise
# set the initial conditions
Tᵢ(z) = Tf - 1e-4 #+ 1e-5 * Ξₜ(z)

set!(model, u=0, v=0, w=0, T=Tᵢ, n1 = n, n2 = n, n3 = n, n4 = n, n5 = n, n6 = n, n7 = n, n8 = n, n9 = n, n10 = n, S=34.5)

simulation = Simulation(model, Δt=1.0, stop_time=1hours)
# Create the simulation

# Add the callback to enforce non-negativity
#add_callback!(simulation, enforce_nonnegative_tracer, IterationInterval(1))
add_callback!(simulation, progress, IterationInterval(1))

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