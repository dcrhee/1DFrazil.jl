using CairoMakie
using JLD2
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat # these constants can be called inside any function
include("test_functions.jl")
# change to number density * crystal volume
# turn off rise

function find_z_mean(var)
    # finds the z-mean of the array
    return vec(mean(var[1, 1, :, :], dims = 1))
end

function find_C_const(indx)
    # finds the z-mean of the array
    return 2*π*Rs[indx]^3/(aspect_ratio*volume)
end

Rs = [0.01, 0.05, 0.15, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 2] .* 1e-3
volume = 1
aspect_ratio = 50

end_name = "match_Feltham_no_collisions"

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

ΔT = time_series.T[1, 1, :, :] .- time_series.T[1, 1, :, 1]
ΔS =  time_series.S[1, 1, :, :] .- time_series.S[1, 1, :, 1]

tmax = maximum(times)
tmin = 0

fig = Figure(size = (850, 850))

ax_ΔT = Axis(fig[1, 1:2];
            ylabel = "ΔT (ᵒC)",
            xlabel = "time (s)",
            limits = ((tmin, tmax), nothing))

ax_C = Axis(fig[1, 3:4];
            ylabel = "C",
            xlabel = "time (s)",
            limits = ((tmin, tmax), nothing))

ax_T = Axis(fig[2, 1:2];
              ylabel = "T-Tf (ᵒC)",
              xlabel = "time (s)",
              limits = ((tmin, tmax), nothing))

ax_n = Axis(fig[2, 3:4];
              ylabel = "n",
              xlabel = "time (s)",
              yscale = log10,
              limits = ((tmin, tmax), nothing))

lines!(ax_ΔT, times, find_z_mean(time_series.T) .- find_z_mean(time_series.T[:, :, :, 1]))
lines!(ax_T, times, find_z_mean(time_series.T) .- Tf)
lines!(ax_C, times, find_C_const(1)*find_z_mean(time_series.n₁), label = "C₁")
lines!(ax_C, times, find_C_const(2)*find_z_mean(time_series.n₂), label = "C₂")
lines!(ax_C, times, find_C_const(3)*find_z_mean(time_series.n₃), label = "C₃")
lines!(ax_C, times, find_C_const(4)*find_z_mean(time_series.n4), label = "C₄")
lines!(ax_C, times, find_C_const(5)*find_z_mean(time_series.n5), label = "C₅")
lines!(ax_C, times, find_C_const(6)*find_z_mean(time_series.n6), label = "C₆")
lines!(ax_C, times, find_C_const(7)*find_z_mean(time_series.n7), label = "C₇")
lines!(ax_C, times, find_C_const(8)*find_z_mean(time_series.n8), label = "C₈")
lines!(ax_C, times, find_C_const(9)*find_z_mean(time_series.n9), label = "C₉")
lines!(ax_C, times, find_C_const(10)*find_z_mean(time_series.n10), label = "C₁₀")

lines!(ax_n, times, find_z_mean(time_series.n₁), label = "n₁")
lines!(ax_n, times, find_z_mean(time_series.n₂), label = "n₂")
lines!(ax_n, times, find_z_mean(time_series.n₃), label = "n₃")
lines!(ax_n, times, find_z_mean(time_series.n4), label = "n₄")
lines!(ax_n, times, find_z_mean(time_series.n5), label = "n₅")
lines!(ax_n, times, find_z_mean(time_series.n6), label = "n₆")
lines!(ax_n, times, find_z_mean(time_series.n7), label = "n₇")
lines!(ax_n, times, find_z_mean(time_series.n8), label = "n₈")
lines!(ax_n, times, find_z_mean(time_series.n9), label = "n₉")
lines!(ax_n, times, find_z_mean(time_series.n10), label = "n₁₀")


axislegend(ax_C)
axislegend(ax_n)
#axislegend(ax_ΔT)
#axislegend(ax_T)

fig

frames = 1:length(times)

save("match_Feltham_no_collisions.png", fig)    