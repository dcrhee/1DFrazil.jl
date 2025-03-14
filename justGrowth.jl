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
depth = 1
numz = 1 #128
grid = RectilinearGrid(size=1, z = (-1, 0), topology=(Flat, Flat, Periodic))
volume = 1
numSizeClasses = 50

# run with the higher values of epsilon and see if it makes a difference using the different formulations

# constants/parameters

# choose radius intervals
aspect_ratio = 50

#
Rs = range(0.01, 2, numSizeClasses) .* 1e-3

#Rs = [0.01, 0.05, 0.15, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 2] .* 1e-3
Hs = 2*Rs/aspect_ratio
Vs = π * Rs.^2 .* Hs

V₁ = Vs[1]
Vₙ = Vs[end]

Tdiff = 1e-3

ϵ = 1e-2#7.4 * 1e-6 #10^-3 # m²s⁻³
ν = 1.95 * 1e-6
Cᵢₙ =  4 * 1e-8
nmax = 10^20 #6429774.231059961 #10^20#10^3 * volume
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

end_name = "just_growth_N" * string(numSizeClasses) * "_T_" * string(Tdiff)

function find_Vi(R, H)
    return π*R^2*H
end

function find_density(T, S)
    ρ = ρₒ * (1 + 7.86 * 1e-4 * (S - 34.5) - 3.87 * 1e-5 * (T + 2) )
    return ρ
end

function find_growth_rate(T, Rᵢ, H)
    #G = kl*Nu/(ρᵢ*Lat) * (Tf - T) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
    G = kl*Nu/(ρᵢ*Lat) * (Tf - T) * 2π * H
    return G
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

function T_forcing_func(z, t, T, S, tracers...)
    tracers_named = NamedTuple{Tuple(Symbol("n$i") for i in 1:numSizeClasses)}(tracers)
    return sum(temperature_forcing_constant(T, S, i) * tracers_named[Symbol("n$i")] for i in 1:numSizeClasses)
end


function nend_forcing_func(i, j, k, grid, clock, model_fields, indx)
    Vᵢ = find_Vi(Rs[indx], Hs[indx])
    # Compute derivative for nᵢ
    nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
    
    # growth term
    Gₙ₋₁ = find_growth_rate(model_fields.T, Rs[indx-1], Hs[indx-1])
    Gₙ = find_growth_rate(model_fields.T, Rs[indx], Hs[indx])
    Vₙ₋₁ = find_Vi(Rs[indx-1], Hs[indx-1])
    Vᵢ = find_Vi(Rs[indx], Hs[indx])

    if Gₙ₋₁[i, j, k] > 0 # growth
        nₙ₋₁ = getfield(model_fields, Symbol("n$(indx-1)"))
        return @inbounds Gₙ₋₁[i, j, k]*nₙ₋₁[i, j, k]/(Vᵢ - Vₙ₋₁)
    else # melt
        final_conc = (Gₙ[i, j, k]*nᵢ[i, j, k])/(Vᵢ - Vₙ₋₁)
        return @inbounds  final_conc
    end

end

function nintermediate_forcing_func(i, j, k, grid, clock, model_fields, indx)
    #ρ = find_density(model_fields.T, model_fields.S)
    
    # Compute derivative for nᵢ
    nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
    nᵢ = max.(nᵢ, 0)

    # Extend velocity array (avoid index errors)
    
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

    
    
    # Apply logic for growth and melting
    if G_im1[i, j, k] > 0  # Growth case
        return @inbounds - (Gᵢ[i, j, k] * nᵢ[i, j, k] / (V_ip1 - Vᵢ) - G_im1[i, j, k] * n_im1[i, j, k] / (Vᵢ - V_im1)) 
    else  # Melt case
        return @inbounds - (G_ip1[i, j, k] * n_ip1[i, j, k] / (V_ip1 - Vᵢ) - Gᵢ[i, j, k] * nᵢ[i, j, k] / (Vᵢ - V_im1)) 
    end
end


function n1_forcing_func(i, j, k, grid, clock, model_fields, indx)
    nᵢ = getfield(model_fields, Symbol("n$indx"))  # Dynamically get field `nᵢ`
    nᵢ = max.(nᵢ, 0)
    
    # growth term
    G₁ = find_growth_rate(model_fields.T, Rs[indx], Hs[indx])
    G₂ = find_growth_rate(model_fields.T, Rs[indx+1], Hs[indx+1])
    V₁ = find_Vi(Rs[indx], Hs[indx])
    V₂ = find_Vi(Rs[indx+1], Hs[indx+1])
    
    if G₁[i, j, k] > 0 # growth
        return @inbounds  - G₁[i, j, k]*model_fields.n1[i, j, k]/(V₂ - V₁)
    else # melt
        return @inbounds   - (G₂[i, j, k]*model_fields.n2[i, j, k]/(V₂ - V₁) - G₁[i, j, k]*model_fields.n1[i, j, k]/V₁) 
    end
end

################################ define forcing functions ###########################################

# Define the range of tracers (e.g., n2 to n100)
n_range_forcing = 2:numSizeClasses-1
n_range_tracers = 2:numSizeClasses
# Create the forcing dictionary dynamically
forcing_dict = Dict(Symbol("n$n") => Forcing(nintermediate_forcing_func, discrete_form=true, parameters=n) for n in n_range_forcing)
n1_forcing = Forcing(n1_forcing_func, discrete_form=true, parameters = 1)
nend_forcing = Forcing(nend_forcing_func, discrete_form=true, parameters = numSizeClasses)

T_field_dependencies = (:T, :S, (Symbol("n$i") for i in 1:numSizeClasses)...)
T_forcing = Forcing(T_forcing_func, field_dependencies=T_field_dependencies)
#S_forcing = Forcing(S_forcing_func, field_dependencies=(:T, :S, :n1, :n2, :n3))

# Merge `n1_forcing` with the dynamic forcing dictionary and convert to NamedTuple
forcing_combined = (; Dict(:T => T_forcing)..., Dict(:n1 => n1_forcing)..., forcing_dict..., Dict(Symbol("n$numSizeClasses") => nend_forcing)...)

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
Tᵢ = Tf - Tdiff#* Ξₜ(z)
#set!(model, u=0, v=0, w=0, T=Tᵢ, n1 = nᵢₙ[1], n2 = nᵢₙ[2], n3 = nᵢₙ[3], n4 = nᵢₙ[4], n5 = nᵢₙ[5], n6 = nᵢₙ[6], n7 = nᵢₙ[7], n8 = nᵢₙ[8], n9 = nᵢₙ[9], n10 = nᵢₙ[10], S=34.5)
# Set initial conditions using `set!`
set!(model, ; u=0, v=0, w=0, T=Tᵢ, S=34.5, ninitial...)

########################## run model #######################################

simulation = Simulation(model, Δt=1.0, stop_time=3hours)

# Define the enforce_nonnegative_tracer function
function enforce_nonnegative_tracer(simulation)
    for n in 1:numSizeClasses
        tracer_data = getproperty(simulation.model.tracers, Symbol("n$n")).data
        @inbounds tracer_data .= max.(tracer_data, 0)
    end
    return nothing
end

# Create the simulation
grid = RectilinearGrid(size=(32, 32, 32), extent=(1, 1, 1))

# Add the callback to enforce non-negativity
add_callback!(simulation, enforce_nonnegative_tracer, IterationInterval(1))
add_callback!(simulation, progress, IterationInterval(20))

conjure_time_step_wizard!(simulation, cfl=1.0, max_Δt=0.01minute)#0.01minute)

output_interval = 0.1minutes

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
#endh