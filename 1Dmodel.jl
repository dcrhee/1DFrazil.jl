using Oceananigans
using CairoMakie
using Oceananigans.Units: minute, minutes, hours
using Oceananigans.BuoyancyModels: g_Earth

using Oceananigans.AbstractOperations: ∂z
using Printf
using Statistics
using Oceananigans.AbstractOperations
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat

# need to fix the forcing functions to include temperature dependence on density
# defauly boundary conditions are no flux

# setup grid: choose 128 data points
grid = RectilinearGrid(size=128, z=(-20, 0), topology=(Flat, Flat, Bounded))

R₁ = 0.001
R₂ = 0.002

Uₐ =  5 # m/s specify the wind strength
Fetch = 1500 # m wind fetch
Tₐ = -20 # atmosphere temperature

end_name = "_Ta_" * string(Tₐ) * "_X_" * string(Fetch) * "_Ua_" * string(Uₐ)

Hₛ = 0.256*Uₐ^2/g_Earth*(1 - (1 + 0.006*(g_Earth*Fetch/Uₐ^2)^(1/2))^(-2))
ω = 2π/24.935*g_Earth/Uₐ*(1.6*Uₐ^2/(g_Earth*Hₛ))^0.625

# specify the waves
wavenumber = ω^2/g_Earth
amplitude = Hₛ/2 # m
wavelength = 2π / wavenumber  # m
frequency = sqrt(g_Earth * wavenumber) # s⁻¹

# The vertical scale over which the Stokes drift of a monochromatic surface wave
# decays away from the surface is `1/2wavenumber`, or
vertical_scale = wavelength / 4π

# Stokes drift velocity at the surface
Uˢ = amplitude^2 * wavenumber * frequency # m s⁻¹
uˢ(z) = Uˢ * exp(z / vertical_scale)
∂z_uˢ(z, t) = 1 / vertical_scale * Uˢ * exp(z / vertical_scale)

# specify the boundary conditions (using the same approach as Near-Inertial Waves and Turbulence Driven by the Growth of Swell, Wagner)
Qᵘ = -ρₐ/ρₒ*Cd*Uₐ^2 # m² s⁻², surface kinematic momentum flux (in time with https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/langmuir-turbulence-in-the-ocean/638FD0E368140E5972144348DB930A38)
u_boundary_conditions = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))
n₁_boundary_conditions = FieldBoundaryConditions(top = GradientBoundaryCondition(0))

Qʰ = 40.0*(Tf - Tₐ)  # W m⁻², surface _heat_ flux
Qᵀ = Qʰ / (ρₒ * cᴾ) # K m s⁻¹, surface _temperature_ flux

u_rise_frazil = 0.001 #sqrt(2*g_Earth*frazil_volume/(crystal_drag*crystal_area)*(1-ρᵢ))

#G2 = p.kl*p.Nu/(p.ρᵢ*p.Lat)*(p.Tf - T) * 2π* p.R₂ * 1/(0.9002 - 0.2634*log(p.H/(2*p.R₂)))
#G1 = p.kl*p.Nu/(p.ρᵢ*p.Lat)*(p.Tf - T) * 2π* p.R₁ * 1/(0.9002 - 0.2634*log(p.H/(2*p.R₁)))

#@inline top_flux(z, t, p) = 40.0*(p.Tf - p.Tₐ)/(p.ρₒ * p.cᴾ)
#@inline top_flux(z, t) = 40.0*(-1.6378 - -20)/(1035 * 3991)
#temptopcondition = FluxBoundaryCondition(top_flux, parameters = (Tf = Tf, Tₐ = Tₐ, ρₒ = ρₒ, cᴾ = cᴾ))
#temptopcondition = FluxBoundaryCondition(top_flux)
#temptopcondition = FluxBoundaryCondition(top_flux, field_dependencies=:T, parameters = (Tf = Tf, Tₐ = Tₐ, ρₒ = ρₒ, cᴾ = cᴾ))

#dTdz = 0.01 # K m⁻¹, impose a gradient boundary condition at the bottom of the domain
T_boundary_conditions = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵀ),
#                                bottom = GradientBoundaryCondition(dTdz))
bottom = ValueBoundaryCondition(Tf))

#T_boundary_conditions = FieldBoundaryConditions(
#bottom = ValueBoundaryCondition(Tf))

#T_boundary_conditions = FieldBoundaryConditions(top = temptopcondition,
#                                bottom = GradientBoundaryCondition(dTdz))
#bottom = ValueBoundaryCondition(Tf))

S_boundary_conditions = FieldBoundaryConditions()

coriolis = FPlane(f=-1.4e-4) # s⁻¹

# Forcing that depends on a `NamedTuple` of parameters `p`
T_forcing_func(z, t, T, n₁, p) = p.k * p.Nu / (p.cᴾ*p.ρ)*(p.Tf - T) * 2π* p.Rₚ * 1/(0.9002 - 0.2634*log(p.H/(2*p.Rₚ))) * n₁
S_forcing_func(z, t, S, T, n₁, p) = S * (1-p.α) * p.ρᵢ/p.ρ * p.k * p.Nu / (p.cᴾ*p.ρ)*(p.Tf - T) * 2π* p.Rₚ * 1/(0.9002 - 0.2634*log(p.H/(2*p.Rₚ))) * n₁
#n1_forcing_func(z, t, T, n₁, p) = - p.k*p.Nu/(p.ρᵢ*p.Lat)*(p.Tf - T) * 2π * (p.R₁ * 1/(0.9002 - 0.2634*log(p.H/(2*p.R₁)))*n₁)/(π*p.H*(p.R₁^2)) #∂z(n₁)
#n2_forcing_func(z, t, T, n₁, n₂, p) = - p.k*p.Nu/(p.ρᵢ*p.Lat)*(p.Tf - T) * 2π * (p.R₂ * 1/(0.9002 - 0.2634*log(p.H/(2*p.R₂)))*n₂ - p.R₁ * 1/(0.9002 - 0.2634*log(p.H/(2*p.R₁)))*n₁)/(π*p.H*(p.R₂^2 - p.R₁^2)) #∂z(n₁)
#n3_forcing_func(z, t, T, n₂, p) = p.k*p.Nu/(p.ρᵢ*p.Lat)*(p.Tf - T) * 2π * (p.R₂ * 1/(0.9002 - 0.2634*log(p.H/(2*p.R₂)))*n₂ )/(π*p.H*(p.R₃^2 - p.R₂^2)) #∂z(n₁)

function find_Vi(H, R)
    return π*R^2*H
end

function find_density(T, S, ρ₀, β, α)
    ρ = ρ₀ * (1.0256550500000001 - α*T + β*S)
    return ρ
end

function find_steady_velocity(Cd, Hᵢ, Rᵢ, ρᵢ, T, S, ρ₀, β, α)
    ρ =  find_density(T, S, ρ₀, β, α)
    Sᵢ = π*Rᵢ^2
    Vᵢ = find_Vi(Hᵢ, Rᵢ)
    u₁ = sqrt.( (2*g_Earth*Vᵢ) / (Cd*Sᵢ) * (1 - ρᵢ/ρ))
    return uᵢ
end

function find_growth_rate()
    G = k1*Nu/(ρᵢ*Lat) * (Tf - T) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
    return G
end

function n1_forcing_func(i, j, k, grid, clock, model_fields, u₁)
    dn_dz = -u₁*∂z(model_fields.n₁)

    ρ = find_density(model_fields.T, model_fields.S, 1027, 0.0078, 1.67*10^(-4))

    return @inbounds dn_dz[i, j, k]
end

n1_forcing = Forcing(n1_forcing_func, discrete_form=true, parameters = 0.001)

T_forcing = Forcing(T_forcing_func, parameters=(cᴾ = cᴾ, k = kl, Nu = Nu, ρ=ρₒ, Rₚ = R₁, H = 0.0004, Tf = Tf), field_dependencies=(:T, :n₁))

S_forcing = Forcing(S_forcing_func, parameters=(cᴾ = cᴾ, α = α, k = kl, Nu = Nu, ρ=ρₒ, ρᵢ=ρᵢ, Rₚ = R₁, H = 0.0004, Tf = Tf), field_dependencies=(:S, :T, :n₁))
#C_forcing = Forcing(C_forcing_func,  parameters=(ufrazil = u_rise_frazil), field_dependencies=(:C))
#n1_forcing = Forcing(n1_forcing_func, parameters=(k = kl, Nu = Nu, ρᵢ=ρᵢ, Lat = Lat, H = 0.0004, R₁ = R₁, R₂ = R₂, Tf = Tf), field_dependencies=(:T, :n₁))
#n2_forcing = Forcing(n2_forcing_func, parameters=(k = kl, Nu = Nu, ρᵢ=ρᵢ, Lat = Lat, H = 0.0004, R₁ = R₁, R₂ = R₂, Tf = Tf), field_dependencies=(:T, :n₁, :n₂))
#n3_forcing = Forcing(n3_forcing_func, parameters=(k = kl, Nu = Nu, ρᵢ=ρᵢ, Lat = Lat, H = 0.0004, R₃ = R₁, R₂ = R₂, Tf = Tf), field_dependencies=(:T, :n₂))

model = NonhydrostaticModel(; grid, coriolis,
advection = WENO(),
timestepper = :RungeKutta3,
tracers = (:T, :S, :n₁, :n₂, :n₃),
buoyancy = SeawaterBuoyancy(),
closure = SmagorinskyLilly(Pr = 1, Cb = 1 / 1),
forcing=(n₁=n1_forcing, S=S_forcing, T=T_forcing),
stokes_drift = UniformStokesDrift(∂z_uˢ=∂z_uˢ),
boundary_conditions = (T=T_boundary_conditions, S=S_boundary_conditions, n₁ = n₁_boundary_conditions))


u, v, w = model.velocities

# set the initial conditions
Ξ(z) = randn() * exp(z / 4)

Ξₜ(z) = randn()  # noise
#Tᵢ(z) = Tf + 1e-6 * Ξₜ(z)
Tᵢ(z) = Tf + 0.01 * Ξₜ(z)

u★ = sqrt(abs(Qᵘ))
uᵢ(z) = u★ * 1e-1 * Ξ(z)
wᵢ(z) = u★ * 1e-1 * Ξ(z)
#Cᵢ(z) = abs.(Ξ(z))*100

width = 5
nᵢ(z) = 1*exp(-(z+10)^2 / (2width^2))

#wᵢ = rand(size(w)...)
#wᵢ .-= mean(wᵢ)

set!(model, w=wᵢ, T=Tᵢ, n₁ = nᵢ, n₂ = nᵢ, n₃ = nᵢ, S=35)

#dz_C = ∂z(C)

simulation = Simulation(model, Δt=45.0, stop_time=8hours)

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

output_interval = 5minutes

fields_to_output = merge(model.velocities, model.tracers, (; νₑ=model.diffusivity_fields.νₑ))

simulation.output_writers[:fields] =
    JLD2OutputWriter(model, fields_to_output,
                    schedule = TimeInterval(output_interval),
                    filename = "1D_fields" * end_name * ".jld2",
                    overwrite_existing = true)

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
n₁ = model.tracers.n₁
n₂ = model.tracers.n₂
n₃ = model.tracers.n₃

run!(simulation)