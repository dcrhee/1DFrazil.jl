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

include("SetupFunctions.jl")
using .SetupFunctions

include("riseVelocityFunctions.jl")
using .riseVelocityFunctions

include("collisionVelocity.jl")
using .collisionVelocity

include("setupSimulation.jl")
using .setupSimulation

using MAT
using JLD2

file = matopen("/Users/cotton/Documents/DPhil/Polynas/Code/non_eqm/Rotation/Full 3D/matlab_full_rot/Code/ODE_solve_interpolate/fall_off/arc/collisionEfficiencyScatteredInterpolate.mat")  # open the .mat file
collision_efficiency = read(file, "collisionEfficiencyScatteredInterpolate")  # read a specific variable from the file
close(file)

# setup grid: choose 128 data points
depth = 0.20
numz = 1 #128
grid = RectilinearGrid(size=1, z = (-1, 0), topology=(Flat, Flat, Periodic))
volume = 1
numSizeClasses = 200

# choose radius intervals
aspect_ratio = 50
Rs = range(0.01, 2, numSizeClasses) .* 1e-3

Hs = 2*Rs/aspect_ratio
Vs = π * Rs.^2 .* Hs

#Vs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #* 1e-9
#Rs = ( (Vs * aspect_ratio) / (2 * π) ) .^(1/3) 
#Hs = 2*Rs/aspect_ratio

V₁ = Vs[1]
Vₙ = Vs[end]

ϵ = 1e-2#7.4 * 1e-6 #10^-3 # m²s⁻³
ν = 1.95 * 1e-6
Cᵢₙ =  4 * 1e-8
nmax = 10^20
ζ = 1 # number of new crystals formed per collision

# work out the indices of the class to count crystal collisions from
Vrem = ζ * V₁

αs, βs, αVolconst, βVolconst  = find_β_α_Vol_constants(Vs, Vrem, V₁, ζ)
const effective_areas = find_mean_collision_radius(Rs, aspect_ratio)
#collision_efficiency = find_collision_efficiency(Rs)
#vgs = find_steady_velocity_via_iteration(Rs, ν, ρₒ, ρᵢ, grav, aspect_ratio)

#@save "vgs200.jld2" vgs  # Saves variable `vgs`
@load "vgs200.jld2" vgs  # Loads it back

coriolis = FPlane(f=-1.4e-4) # s⁻¹

cvelnum = 1 # 1 is old cylinder, 2 is new cylinder, 3 is new spherical
pnum = 2 # 1 is mean n, 2 is sum over nj

end_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)

collision_velocity_parameterisation_num = cvelnum # 1 is old cylinder, 2 is new cylinder, 3 is new spherical
concentration_parameterisation_num = pnum # 1 is mean n, 2 is sum over nj
crystal_size_collision_redistribution = 2 # 1 is old redistribution, 2 is new redistribution
effective_radius = true # use their equivalent radius
averaged_radius = false # use their averaged radius
max_radius = false #true

if collision_velocity_parameterisation_num == 2
    end_name = end_name * "_new_cyl"
elseif collision_velocity_parameterisation_num == 3
    end_name = end_name * "_new_spherical"
end
if concentration_parameterisation_num == 2
    end_name = end_name * "_sum_nj"
end
if averaged_radius
    end_name = end_name * "_r_av"
end
if max_radius
    end_name = end_name * "_r_max"
end
if crystal_size_collision_redistribution == 2
    end_name = end_name * "_redistribute"
end

# setup model 

# pre calculate all the turbulent velocities, all the rise velocities and all the relative velocities
include("turbulentVelocityFunctions.jl")
using .turbulentVelocityFunctions

if effective_radius
    if collision_velocity_parameterisation_num == 2 # 2 is new cylinder
        get_vₜ = get_vₜ_effective_radius_collision_velocity_parameterisation_num_2
        get_vₜ_same_radius = get_vₜ_same_radius_effective_radius_collision_velocity_parameterisation_num_2
    else # 1 is old cylinder, 3 is new spherical
        get_vₜ = get_vₜ_effective_radius_collision_velocity_parameterisation_num_13
        get_vₜ_same_radius = get_vₜ_same_radius_effective_radius_collision_velocity_parameterisation_num_13
    end
elseif averaged_radius
    if collision_velocity_parameterisation_num == 2 # 2 is new cylinder
        get_vₜ = get_vₜ_averaged_radius_collision_velocity_parameterisation_num_2
        get_vₜ_same_radius = get_vₜ_same_radius_averaged_radius_collision_velocity_parameterisation_num_2
    else # 1 is old cylinder, 3 is new spherical
        get_vₜ = get_vₜ_averaged_radius_collision_velocity_parameterisation_num_13
        get_vₜ_same_radius = get_vₜ_same_radius_averaged_radius_collision_velocity_parameterisation_num_13
    end
else
    if collision_velocity_parameterisation_num == 2 # 2 is new cylinder
        get_vₜ = get_vₜ_collision_velocity_parameterisation_num_2
        get_vₜ_same_radius = get_vₜ_same_radius_collision_velocity_parameterisation_num_2
    else # 1 is old cylinder, 3 is new spherical
        get_vₜ = get_vₜ_collision_velocity_parameterisation_num_13
        get_vₜ_same_radius = get_vₜ_same_radius_collision_velocity_parameterisation_num_13
    end
end

# get a matrix of the relative velocities
if concentration_parameterisation_num == 1
    const vcoll_matrix = zeros(length(Rs))
    for Riindx in eachindex(Rs)
        vₜ = get_vₜ_same_radius(Riindx, ϵ, ν)
        vᵣ = find_steady_velocity(Riindx)
        vcoll_matrix[Riindx] = find_vcoll(vᵣ, vₜ, collision_velocity_parameterisation_num)
    end
else
    const vcoll_matrix = zeros(length(Rs), length(Rs))
    for Riindx in eachindex(Rs)
        uᵢ = find_steady_velocity(Riindx)
        for Rjindx in eachindex(Rs)
            uⱼ = find_steady_velocity(Rjindx)
            vᵣ = abs(uᵢ - uⱼ)
            vₜ = get_vₜ(Riindx, Rjindx, ϵ, ν)
            vcoll_matrix[Riindx, Rjindx] = find_vcoll(vᵣ, vₜ, collision_velocity_parameterisation_num)
        end
    end
end

# here is fine

include("encounterFrequencyConstants.jl")
using .encounterFrequencyConstants

# choose the right encounter frquency function
if effective_radius
    if collision_velocity_parameterisation_num == 3 # 3 is new spherical
        get_Fenc = get_Fenc_effective_radius_collision_velocity_parameterisation_num_3
        get_Fenc_ndensity = get_Fenc_p1_effective_radius_collision_velocity_parameterisation_num_3
    else # 1 is old cylinder, 2 is new cylinder, 
        get_Fenc = get_Fenc_effective_radius_collision_velocity_parameterisation_num_12
        get_Fenc_ndensity = get_Fenc_p1_effective_radius_collision_velocity_parameterisation_num_12
    end
elseif averaged_radius
    if collision_velocity_parameterisation_num == 3 # 3 is new spherical
        get_Fenc = get_Fenc_average_radius_collision_velocity_parameterisation_num_3
        get_Fenc_ndensity = get_Fenc_p1_average_radius_collision_velocity_parameterisation_num_3
    else # 1 is old cylinder, 2 is new cylinder, 
        get_Fenc = get_Fenc_average_radius_collision_velocity_parameterisation_num_12
        get_Fenc_ndensity = get_Fenc_p1_average_radius_collision_velocity_parameterisation_num_12
    end
else
    if collision_velocity_parameterisation_num == 3 # 3 is new spherical
        get_Fenc =  get_Fenc_collision_velocity_parameterisation_num_3
        get_Fenc_ndensity = get_Fenc_p1_collision_velocity_parameterisation_num_3
    else # 1 is old cylinder, 2 is new cylinder, 
        get_Fenc = get_Fenc_collision_velocity_parameterisation_num_12
        get_Fenc_ndensity = get_Fenc_p1_collision_velocity_parameterisation_num_12
    end
end   

# get a matrix of the relative velocities
if concentration_parameterisation_num == 1
    const encounterFrequencyConstants_matrix = zeros(length(Rs))
    for Riindx in eachindex(Rs)
        encounterFrequencyConstants_matrix[Riindx] = get_Fenc_ndensity(Riindx)
    end
else
    const encounterFrequencyConstants_matrix = zeros(length(Rs), length(Rs))
    for Riindx in eachindex(Rs)
        uᵢ = find_steady_velocity(Riindx)
        for Rjindx in eachindex(Rs)
            encounterFrequencyConstants_matrix[Riindx, Rjindx] = get_Fenc(Riindx, Rjindx)
        end
    end
end

include("collisionFrequencyFunctions.jl")
using .collisionFrequencyFunctions

if concentration_parameterisation_num == 1
    coll_freq_concentration_parameterisation = coll_freq_concentration_parameterisation_num_1
else
    coll_freq_concentration_parameterisation = coll_freq_concentration_parameterisation_num_2
end

################################ define forcing functions ###########################################
include("forcingFunctions.jl")
using .forcingFunctions
#include("forcingFunctionTest.jl")
#using .forcingFunctionTest

# Define the range of tracers (e.g., n2 to n100)
n_range_forcing = 2:numSizeClasses-1
n_range_tracers = 2:numSizeClasses
# Create the forcing dictionary dynamically
forcing_dict = Dict(Symbol("n$n") => Forcing(nintermediate_forcing_func, discrete_form=true, parameters=n) for n in n_range_forcing)
n1_forcing = Forcing(n1_forcing_func, discrete_form=true, parameters = 1)
nend_forcing = Forcing(nend_forcing_func, discrete_form=true, parameters = numSizeClasses)
#T_forcing = Forcing(T_forcing_func, field_dependencies=(:T, :S, :n1, :n2, :n3, :n4, :n5, :n6, :n7, :n8, :n9, :n10))
#S_forcing = Forcing(S_forcing_func, field_dependencies=(:T, :S, :n1, :n2, :n3))

# Merge `n1_forcing` with the dynamic forcing dictionary and convert to NamedTuple
forcing_combined = (; Dict(:n1 => n1_forcing)..., forcing_dict..., Dict(Symbol("n$numSizeClasses") => nend_forcing)...)

######################### setup model ########################################

# Define the model with dynamically generated tracers and forcing
model = NonhydrostaticModel(; grid,
    advection = Centered(), # WENO(),
    timestepper = :RungeKutta3,
    buoyancy = SeawaterBuoyancy(),
    tracers = (:T, :S, :n1, Symbol.(["n$n" for n in n_range_tracers])...),  # Dynamically create the list
    forcing = forcing_combined  # Merge manually defined n1 with the dictionary
)

u, v, w = model.velocities
#nᵢₙ = Cᵢₙ * aspect_ratio * volume / (2 * π * length(Rs)) * 1 ./ Rs.^3
nᵢₙ = 6429774.231059961/length(Rs) * ones(length(Rs))

# Create a dictionary for dynamically setting `n` values
ninitial_dict = Dict(Symbol("n$n") => nᵢₙ[n] for n in 1:numSizeClasses)
ninitial = (; ninitial_dict...)

# set the initial conditions
Tᵢ = Tf - 1e-4 #* Ξₜ(z)
#set!(model, u=0, v=0, w=0, T=Tᵢ, n1 = nᵢₙ[1], n2 = nᵢₙ[2], n3 = nᵢₙ[3], n4 = nᵢₙ[4], n5 = nᵢₙ[5], n6 = nᵢₙ[6], n7 = nᵢₙ[7], n8 = nᵢₙ[8], n9 = nᵢₙ[9], n10 = nᵢₙ[10], S=34.5)
# Set initial conditions using `set!`
set!(model, ; u=0, v=0, w=0, T=Tᵢ, S=34.5, ninitial...)

########################## run model #######################################

#simulation = Simulation(model, Δt=1.0, stop_time=8000minutes)
simulation = Simulation(model, Δt=1.0, stop_time=30minutes)
#simulation = Simulation(model, Δt=1.0, stop_time=10minutes)
#simulation = Simulation(model, Δt=1.0, stop_time=4000minutes) not long enough for epsilon 1e-2 sum nj

# Create the simulation
grid = RectilinearGrid(size=(32, 32, 32), extent=(1, 1, 1))

# Add the callback to enforce non-negativity
add_callback!(simulation, enforce_nonnegative_tracer, IterationInterval(1))
add_callback!(simulation, progress, IterationInterval(1))

#conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=0.000001minute) # for n = 1, p = 1, eps = 1e-2
#conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=0.00001minute) # for n = 1, p = 2, eps = 1e-2

conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=0.001minute) # for n = 1, p = 2
#conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=0.1minute) # for n = 1, p = 2

#conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=0.001minute) # for n = 1, p = 2 for n sum

#conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=0.02minute) # for n = 1, p = 2

#simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))

output_interval = 0.02minutes

#output_interval = 0.1minutes

fields_to_output = merge(model.velocities, model.tracers)

simulation.output_writers[:fields] =
    JLD2OutputWriter(model, fields_to_output,
                    schedule = TimeInterval(output_interval),
                    filename = "1D_fields" * end_name * ".jld2",
                    overwrite_existing = true)

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S

run!(simulation)
#end
#endPh
