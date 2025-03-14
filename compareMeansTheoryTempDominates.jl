using CairoMakie
using JLD2
using PerceptualColourMaps
using Oceananigans
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat # these constants can be called inside any function
using Statistics
include("test_functions.jl")
# change to number density * crystal volume
# turn off rise

function temperature_forcing_constant(T, S, indx)
    # constant in front of each concentration
    #Rᵢ = Rs[indx]
    H = Hs[indx]
    ρ = find_density(T, S)
    Tconstᵢ = kl*Nu/(volume * ρ*cᴾ) * 2π * H
    #Tconstᵢ = kl*Nu/(volume * ρ*cᴾ) * (Tf - T) *2π * Rᵢ *  1/(0.9002 - 0.2634*log(H/(2*Rᵢ)))
    return Tconstᵢ
end

function find_temp_time_const(T, S, n₁, n₂, n₃, n₄, n₅, n₆, n₇, n₈, n₉, n10)
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

function temp_theory(t, Tₒ, S₀, n₁ₒ, n₂ₒ, n₃ₒ, n₄₀, n₅₀, n₆₀, n₇₀, n₈₀, n₉₀, n₁₀₀)
    # estimate for the temperature exponential
    A =  find_temp_time_const(Tₒ, S₀, n₁ₒ, n₂ₒ, n₃ₒ, n₄₀, n₅₀, n₆₀, n₇₀, n₈₀, n₉₀, n₁₀₀)
    T = Tf .- (Tf - Tₒ)*exp.(-A*t)
    return T
end

function n1_theory(t, Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ)
    A =  T_const(Tₒ, Sₒ, n₁ₒ, n₂ₒ, n₃ₒ, R₁, R₂, R₃) # input the concentrations at t = 0
    G₁ = G_const(R₁)
    V₁ = find_Vi(R₁)
    V₂ = find_Vi(R₂)
    n₁ = n₁ₒ*exp.( (1 .- exp.(-A*t)) * G₁/(A*(V₂-V₁)) * (Tₒ - Tf) )
    return n₁
end

function find_z_mean(var)
    # finds the z-mean of the array
    return vec(mean(var[1, 1, :, :], dims = 1))
end

function find_C_const(indx)
    # finds the z-mean of the array
    return 2*π*Rs[indx]^3/(aspect_ratio*volume)
end

aspect_ratio = 50
Rs = [0.01, 0.05, 0.15, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 2] .* 1e-3
Hs = 2*Rs/aspect_ratio
volume = 1

end_name = "match_Feltham_no_collisions"

#colours = cmap("Gouldian", N = length(Rs))


time_series = (;
    w = FieldTimeSeries("1D_fields"* end_name * ".jld2", "w"),
    u = FieldTimeSeries("1D_fields"* end_name * ".jld2", "u"),
    v = FieldTimeSeries("1D_fields"* end_name * ".jld2", "v"),
    T = FieldTimeSeries("1D_fields"* end_name * ".jld2", "T"),
    S = FieldTimeSeries("1D_fields"* end_name * ".jld2", "S"),
    n₁ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n1"),
    n₂ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n2"),
    n₃ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n3"),
    n4 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n4"),
    n5 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n5"),
    n6 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n6"),
    n7 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n7"),
    n8 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n8"),
    n9 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n9"),
    n10 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n10"),
    )

times = time_series.w.times
times = times/(3600*24)


n₁₀ =  time_series.n₁[:, :, :, 1][1]
n₂₀ =  time_series.n₂[:, :, :, 1][1]
n₃₀ =  time_series.n₃[:, :, :, 1][1]
n₄₀ =  time_series.n4[:, :, :, 1][1]
n₅₀ =  time_series.n5[:, :, :, 1][1]
n₆₀ =  time_series.n6[:, :, :, 1][1]
n₇₀ =  time_series.n7[:, :, :, 1][1]
n₈₀ =  time_series.n8[:, :, :, 1][1]
n₉₀ =  time_series.n9[:, :, :, 1][1]
n₁₀₀ =  time_series.n10[:, :, :, 1][1]
Tₒ = find_z_mean(time_series.T[:, :, :, 1])[1]
S₀ = find_z_mean(time_series.S[:, :, :, 1])[1]

Tguess = temp_theory(times * 3600 * 24, Tₒ, S₀, n₁₀, n₂₀, n₃₀, n₄₀, n₅₀, n₆₀, n₇₀, n₈₀, n₉₀, n₁₀₀)

ΔT = time_series.T[1, 1, :, :] .- time_series.T[1, 1, :, 1]

tmax = maximum(times)
tmin = 0

fig = Figure(size = (850, 850))

ax_ΔT = Axis(fig[1, 1:2];
            ylabel = "ΔT (ᵒC)",
            xlabel = "time (days)",
            limits = ((tmin, tmax), nothing))

ax_C = Axis(fig[1, 3:4];
            ylabel = "C",
            xlabel = "time (days)",
            yscale = log10,
            limits = ((tmin, tmax), nothing))

ax_T = Axis(fig[2, 1:2];
            ylabel = "T-Tf (ᵒC)",
            xlabel = "time (days)",
            limits = ((tmin, tmax), nothing))

ax_n = Axis(fig[2, 3:4];
            ylabel = "n",
            xlabel = "time (days)",
            yscale = log10,
            limits = ((tmin, tmax), (0.1, nothing)))

lines!(ax_ΔT, times, find_z_mean(time_series.T) .- find_z_mean(time_series.T[:, :, :, 1]), color=:black, label = "model")
lines!(ax_T, times, find_z_mean(time_series.T) .- Tf, color=:black, label = "model")

lines!(ax_ΔT, times, Tguess .- find_z_mean(time_series.T[:, :, :, 1]), label = "theory",  linestyle = :dash)
lines!(ax_T, times, Tguess .- Tf, label = "theory", linestyle = :dash)

lines!(ax_C, times, find_C_const(1)*find_z_mean(time_series.n₁), label = "C₁", color = colours[1])
lines!(ax_C, times, find_C_const(2)*find_z_mean(time_series.n₂), label = "C₂", color = colours[2])
lines!(ax_C, times, find_C_const(3)*find_z_mean(time_series.n₃), label = "C₃", color = colours[3])
lines!(ax_C, times, find_C_const(4)*find_z_mean(time_series.n4), label = "C₄", color = colours[4])
lines!(ax_C, times, find_C_const(5)*find_z_mean(time_series.n5), label = "C₅", color = colours[5])
lines!(ax_C, times, find_C_const(6)*find_z_mean(time_series.n6), label = "C₆", color = colours[6])
lines!(ax_C, times, find_C_const(7)*find_z_mean(time_series.n7), label = "C₇", color = colours[7])
lines!(ax_C, times, find_C_const(8)*find_z_mean(time_series.n8), label = "C₈", color = colours[8])
lines!(ax_C, times, find_C_const(9)*find_z_mean(time_series.n9), label = "C₉", color = colours[9])
lines!(ax_C, times, find_C_const(10)*find_z_mean(time_series.n10), label = "C₁₀", color = colours[10])

lines!(ax_n, times, find_z_mean(time_series.n₁), label = "n₁", color = colours[1])
lines!(ax_n, times, find_z_mean(time_series.n₂), label = "n₂", color = colours[2])
lines!(ax_n, times, find_z_mean(time_series.n₃), label = "n₃", color = colours[3])
lines!(ax_n, times, find_z_mean(time_series.n4), label = "n₄", color = colours[4])
lines!(ax_n, times, find_z_mean(time_series.n5), label = "n₅", color = colours[5])
lines!(ax_n, times, find_z_mean(time_series.n6), label = "n₆", color = colours[6])
lines!(ax_n, times, find_z_mean(time_series.n7), label = "n₇", color = colours[7])
lines!(ax_n, times, find_z_mean(time_series.n8), label = "n₈", color = colours[8])
lines!(ax_n, times, find_z_mean(time_series.n9), label = "n₉", color = colours[9])
lines!(ax_n, times, find_z_mean(time_series.n10), label = "n₁₀", color = colours[10])


axislegend(ax_C)
axislegend(ax_n)
axislegend(ax_ΔT, position=(:right, :bottom, ))
#axislegend(ax_T)

fig

frames = 1:length(times)

save(end_name * "_compareMeansTheoryTempDominates.png", fig)