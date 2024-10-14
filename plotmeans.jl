using CairoMakie
using JLD2

# change to number density * crystal volume
# turn off rise

Uₐ =  5 # m/s specify the wind strength
Fetch = 1500 # m wind fetch
Tₐ = -20 # atmosphere temperature

end_name = "_Ta_" * string(Tₐ) * "_X_" * string(Fetch) * "_Ua_" * string(Uₐ)

time_series = (;
     w = FieldTimeSeries("1D_fields"* end_name * ".jld2", "w"),
     u = FieldTimeSeries("1D_fields"* end_name * ".jld2", "u"),
     v = FieldTimeSeries("1D_fields"* end_name * ".jld2", "v"),
     T = FieldTimeSeries("1D_fields"* end_name * ".jld2", "T"),
     S = FieldTimeSeries("1D_fields"* end_name * ".jld2", "S"),
     n₁ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n₁"),
     n₂ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n₂"),
     n₃ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n₃"),
     )

times = time_series.w.times

function find_z_mean(var)
    # finds the z-mean of the array
    return vec(mean(var[1, 1, :, :], dims = 1))
end

function find_density(T, S, ρ₀ = 1027.0, β = 0.0078, α = 1.67*10^(-4))
    ρ = ρ₀ * (1.0256550500000001 - α*T + β*S)
    return ρ
end

ΔT = time_series.T[1, 1, :, :] .- time_series.T[1, 1, :, 1]
ΔS =  time_series.S[1, 1, :, :] .- time_series.S[1, 1, :, 1]


fig = Figure(size = (850, 850))

ax_ΔT = Axis(fig[1, 1:2];
            ylabel = "ΔT (ᵒC)",
            xlabel = "time (s)")

ax_ΔS = Axis(fig[1, 3:4];
            ylabel = "ΔS (ppt)",
            xlabel = "time (s)")

ax_rho = Axis(fig[1, 5:6];
            ylabel = "ρ (kg/m³)",
            xlabel = "time (s)")

ax_T = Axis(fig[2, 1:2];
              ylabel = "T (ᵒC)",
              xlabel = "time (s)")

ax_S = Axis(fig[2, 3:4];
              ylabel = "S (ppt)",
              xlabel = "time (s)")
              #limits = ((minimum(time_series.S), maximum(time_series.S)), nothing))

ax_C = Axis(fig[2, 5:6];
              ylabel = "n",
              xlabel = "time (s)")

lines!(ax_ΔT, times, find_z_mean(time_series.T) .- find_z_mean(time_series.T[:, :, :, 1]))
lines!(ax_ΔS, times, find_z_mean(time_series.S) .- find_z_mean(time_series.S[:, :, :, 1]))
lines!(ax_rho, times, find_density(find_z_mean(time_series.T), find_z_mean(time_series.S)))
lines!(ax_T, times, find_z_mean(time_series.T))
lines!(ax_S, times, find_z_mean(time_series.S))
lines!(ax_C, times, find_z_mean(time_series.n₁), label = "n₁")
lines!(ax_C, times, find_z_mean(time_series.n₂), label = "n₂")
lines!(ax_C, times, find_z_mean(time_series.n₃), label = "n₃")

axislegend(ax_C)

fig

frames = 1:length(times)

save("no_rise.png", fig)