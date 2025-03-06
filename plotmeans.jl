using CairoMakie
using JLD2
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat # these constants can be called inside any function
include("test_functions.jl")
# change to number density * crystal volume
# turn off rise

Uₐ =  5 # m/s specify the wind strength
Fetch = 1500 # m wind fetch
Tₐ = -20 # atmosphere temperature

R₁ = 0.001
R₂ = 0.002
R₃ = 0.003

#end_name = "_Ta_" * string(Tₐ) * "_X_" * string(Fetch) * "_Ua_" * string(Uₐ)
end_name = "match_Feltham"


#time_series = (;
#     w = FieldTimeSeries("1D_fields"* end_name * ".jld2", "w"),
#     u = FieldTimeSeries("1D_fields"* end_name * ".jld2", "u"),
#     v = FieldTimeSeries("1D_fields"* end_name * ".jld2", "v"),
#     T = FieldTimeSeries("1D_fields"* end_name * ".jld2", "T"),
#     S = FieldTimeSeries("1D_fields"* end_name * ".jld2", "S"),
#     n₁ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n₁"),
#     n₂ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n₂"),
#     n₃ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n₃"),
#     )

time_series = (;
     w = FieldTimeSeries("1D_fields"* end_name * ".jld2", "w"),
     u = FieldTimeSeries("1D_fields"* end_name * ".jld2", "u"),
     v = FieldTimeSeries("1D_fields"* end_name * ".jld2", "v"),
     T = FieldTimeSeries("1D_fields"* end_name * ".jld2", "T"),
     S = FieldTimeSeries("1D_fields"* end_name * ".jld2", "S"),
     n₁ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n1"),
     n₂ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n2"),
     n₃ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n3"),
     )

times = time_series.w.times

function find_z_mean(var)
    # finds the z-mean of the array
    return vec(mean(var[1, 1, :, :], dims = 1))
end

function find_density(T, S, ρ₀ = 1027.0, β = 0.0078, α = 1.67*10^(-4))
    ρ = ρ₀ * (1.0256550500000001 .- α*T .+ β*S)
    return ρ
end

function temp_theory(t, Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ)
    # estimate for the temperature exponential
    A =  T_const(Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ, R₁, R₂, R₃) # input the concentrations at t = 0
    T = Tf .- (Tf - Tₒ)*exp.(-A*t)
    return T
end

function salinity_theory(t, Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ)
    A =  T_const(Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ, R₁, R₂, R₃) # input the concentrations at t = 0
    S = Sₒ*exp.(S_const(Tₒ)* (exp.(-A*t) .- 1)   )
    return S   
end

function n1_theory(t, Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ)
    A =  T_const(Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ, R₁, R₂, R₃) # input the concentrations at t = 0
    G₁ = G_const(R₁)
    V₁ = find_Vi(R₁)
    V₂ = find_Vi(R₂)
    n₁ = n₁ₒ*exp.( (1 .- exp.(-A*t)) * G₁/(A*(V₂-V₁)) * (Tₒ - Tf) )
    return n₁
end

function n2_theory(t, Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ)
    τ =  1/T_const(Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ, R₁, R₂, R₃) # input the concentrations at t = 0
    G₁ = G_const(R₁)
    G₂ = G_const(R₂)
    V₁ = find_Vi(R₁)
    V₂ = find_Vi(R₂)
    V₃ = find_Vi(R₃)
    
    dV₃ = V₃ - V₂
    dV₂ = V₂ - V₁

    first_term = -dV₃*G₁*n₁ₒ*exp.( (Tₒ - Tf) * τ  .* ( G₁*dV₃*(1 .- exp.(-t/τ)) .+ dV₂ * G₂ * exp.(-t/τ)  ) / (dV₃*dV₂))
    second_term = (-dV₂*G₂*n₂ₒ + dV₃*G₁*(n₁ₒ + n₂ₒ)) * exp(G₂ * (Tₒ - Tf) * τ / dV₃)
    n₂ = 1/(dV₃*G₁ - dV₂*G₂)*exp.(G₂*τ*(Tf - Tₒ)/dV₃ * exp.(-t/τ)) .* (first_term .+ second_term)

    return n₂
end

function n3_theory(t, Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ)
    τ =  1/T_const(Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ, R₁, R₂, R₃) # input the concentrations at t = 0
    G₁ = G_const(R₁)
    G₂ = G_const(R₂)
    V₁ = find_Vi(R₁)
    V₂ = find_Vi(R₂)
    V₃ = find_Vi(R₃)
    dV₃ = V₃ - V₂
    dV₂ = V₂ - V₁

    first_term = dV₂*G₂*n₁ₒ * exp.((1 .- exp.(-t/τ))*G₁*(Tₒ - Tf) *τ/ dV₂ )
    second_term = (dV₂*G₂*n₂ₒ - dV₃*G₁*(n₁ₒ + n₂ₒ)) * exp.((1 .- exp.(-t/τ))*G₂*(Tₒ - Tf) *τ/ dV₃ )
    third_term = (dV₃*G₁ - dV₂*G₂) * (n₁ₒ + n₂ₒ + n₃ₒ)
    n₃ = 1/(dV₃*G₁ - dV₂*G₂) * (first_term .+ second_term .+ third_term)

    return n₃
end

ΔT = time_series.T[1, 1, :, :] .- time_series.T[1, 1, :, 1]
ΔS =  time_series.S[1, 1, :, :] .- time_series.S[1, 1, :, 1]

tmax = 300
tmin = 0

fig = Figure(size = (850, 850))

ax_ΔT = Axis(fig[1, 1:2];
            ylabel = "ΔT (ᵒC)",
            xlabel = "time (s)",
            limits = ((tmin, tmax), nothing))

ax_ΔS = Axis(fig[1, 3:4];
            ylabel = "ΔS (ppt)",
            xlabel = "time (s)",
            limits = ((tmin, tmax), nothing))

ax_rho = Axis(fig[1, 5:6];
            ylabel = "ρ (kg/m³)",
            xlabel = "time (s)",
            limits = ((tmin, tmax), nothing))

ax_T = Axis(fig[2, 1:2];
              ylabel = "T-Tf (ᵒC)",
              xlabel = "time (s)",
              limits = ((tmin, tmax), nothing))

ax_S = Axis(fig[2, 3:4];
              ylabel = "S (ppt)",
              xlabel = "time (s)",
              limits = ((tmin, tmax), nothing))

ax_C = Axis(fig[2, 5:6];
              ylabel = "n",
              xlabel = "time (s)",
              limits = ((tmin, tmax), nothing))

Tguess = temp_theory(times, find_z_mean(time_series.T[:, :, :, 1])[1], find_z_mean(time_series.S[:, :, :, 1])[1], find_z_mean(time_series.n₁[:, :, :, 1])[1], find_z_mean(time_series.n₂[:, :, :, 1])[1], find_z_mean(time_series.n₃[:, :, :, 1])[1])
Sguess = salinity_theory(times, find_z_mean(time_series.T[:, :, :, 1])[1], find_z_mean(time_series.S[:, :, :, 1])[1], find_z_mean(time_series.n₁[:, :, :, 1])[1], find_z_mean(time_series.n₂[:, :, :, 1])[1], find_z_mean(time_series.n₃[:, :, :, 1])[1])
n₁guess = n1_theory(times, find_z_mean(time_series.T[:, :, :, 1])[1], find_z_mean(time_series.S[:, :, :, 1])[1], find_z_mean(time_series.n₁[:, :, :, 1])[1], find_z_mean(time_series.n₂[:, :, :, 1])[1], find_z_mean(time_series.n₃[:, :, :, 1])[1])
n₂guess = n2_theory(times, find_z_mean(time_series.T[:, :, :, 1])[1], find_z_mean(time_series.S[:, :, :, 1])[1], find_z_mean(time_series.n₁[:, :, :, 1])[1], find_z_mean(time_series.n₂[:, :, :, 1])[1], find_z_mean(time_series.n₃[:, :, :, 1])[1])
n₃guess = n3_theory(times, find_z_mean(time_series.T[:, :, :, 1])[1], find_z_mean(time_series.S[:, :, :, 1])[1], find_z_mean(time_series.n₁[:, :, :, 1])[1], find_z_mean(time_series.n₂[:, :, :, 1])[1], find_z_mean(time_series.n₃[:, :, :, 1])[1])

lines!(ax_ΔT, times, find_z_mean(time_series.T) .- find_z_mean(time_series.T[:, :, :, 1]))
lines!(ax_ΔS, times, find_z_mean(time_series.S) .- find_z_mean(time_series.S[:, :, :, 1]))
lines!(ax_rho, times, find_density(find_z_mean(time_series.T), find_z_mean(time_series.S)))
lines!(ax_T, times, find_z_mean(time_series.T) .- Tf)
lines!(ax_S, times, find_z_mean(time_series.S))
lines!(ax_C, times, find_z_mean(time_series.n₁), label = "n₁")
lines!(ax_C, times, find_z_mean(time_series.n₂), label = "n₂")
lines!(ax_C, times, find_z_mean(time_series.n₃), label = "n₃")
#lines!(ax_C, times, find_z_mean(time_series.n₁) .+ find_z_mean(time_series.n₂) .+ find_z_mean(time_series.n₃), label = "nₜ")

# add in guesses
lines!(ax_ΔT, times, Tguess .- find_z_mean(time_series.T[:, :, :, 1]), label = "theory",  linestyle = :dash)
lines!(ax_T, times, Tguess .- Tf, label = "theory", linestyle = :dash)
lines!(ax_ΔS, times, Sguess .- find_z_mean(time_series.S[:, :, :, 1]), label = "theory",  linestyle = :dash)
lines!(ax_S, times, Sguess, label = "theory", linestyle = :dash)
lines!(ax_C, times, n₁guess, label = "n₁ theory", linestyle = :dash)
lines!(ax_C, times, n₂guess, label = "n₂ theory", linestyle = :dash)
lines!(ax_C, times, n₃guess, label = "n₃ theory", linestyle = :dash)

axislegend(ax_C)
axislegend(ax_ΔT)
axislegend(ax_T)

fig

frames = 1:length(times)

save("no_coll_mixed.png", fig)    