using CairoMakie
using JLD2
using PerceptualColourMaps
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

feltham_name = "match_Feltham"


#colours = cmap("Gouldian", N = length(Rs))
colours = cmap("CBTL1", N = Int(1.2*length(Rs))) 
end_name = "match_Feltham_no_n_max"
new_label = "n sum"
new_label = "new_spherical"
new_label = "no nmax"

feltham_time_series = (;
     w = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "w"),
     u = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "u"),
     v = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "v"),
     T = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "T"),
     S = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "S"),
     n₁ = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n1"),
     n₂ = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n2"),
     n₃ = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n3"),
     n4 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n4"),
     n5 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n5"),
     n6 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n6"),
     n7 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n7"),
     n8 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n8"),
     n9 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n9"),
     n10 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n10"),
     )

feltham_times = feltham_time_series.w.times
feltham_times = feltham_times/(3600*24)

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

ΔT = time_series.T[1, 1, :, :] .- time_series.T[1, 1, :, 1]
ΔS =  time_series.S[1, 1, :, :] .- time_series.S[1, 1, :, 1]

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
              limits = ((tmin, tmax), nothing))

lines!(ax_ΔT, feltham_times, find_z_mean(feltham_time_series.T) .- find_z_mean(time_series.T[:, :, :, 1]), linestyle=:dash, color=:black, label = "H&F")
lines!(ax_T, feltham_times, find_z_mean(feltham_time_series.T) .- Tf, linestyle=:dash, color=:black, label = "H&F")
lines!(ax_ΔT, times, find_z_mean(time_series.T) .- find_z_mean(time_series.T[:, :, :, 1]), color=:black, label = new_label)
lines!(ax_T, times, find_z_mean(time_series.T) .- Tf, color=:black, label = new_label)

lines!(ax_C, feltham_times, find_C_const(1)*find_z_mean(feltham_time_series.n₁), linestyle=:dash,  color = colours[1])
lines!(ax_C, feltham_times, find_C_const(2)*find_z_mean(feltham_time_series.n₂), linestyle=:dash,  color = colours[2])
lines!(ax_C, feltham_times, find_C_const(3)*find_z_mean(feltham_time_series.n₃), linestyle=:dash,  color = colours[3])
lines!(ax_C, feltham_times, find_C_const(4)*find_z_mean(feltham_time_series.n4), linestyle=:dash,  color = colours[4])
lines!(ax_C, feltham_times, find_C_const(5)*find_z_mean(feltham_time_series.n5), linestyle=:dash,  color = colours[5])
lines!(ax_C, feltham_times, find_C_const(6)*find_z_mean(feltham_time_series.n6), linestyle=:dash,  color = colours[6])
lines!(ax_C, feltham_times, find_C_const(7)*find_z_mean(feltham_time_series.n7), linestyle=:dash,  color = colours[7])
lines!(ax_C, feltham_times, find_C_const(8)*find_z_mean(feltham_time_series.n8), linestyle=:dash,  color = colours[8])
lines!(ax_C, feltham_times, find_C_const(9)*find_z_mean(feltham_time_series.n9), linestyle=:dash,  color = colours[9])
lines!(ax_C, feltham_times, find_C_const(10)*find_z_mean(feltham_time_series.n10), linestyle=:dash, color = colours[10])

lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n₁), linestyle=:dash, color = colours[1])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n₂), linestyle=:dash,  color = colours[2])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n₃), linestyle=:dash,  color = colours[3])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n4), linestyle=:dash,  color = colours[4])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n5), linestyle=:dash,  color = colours[5])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n6), linestyle=:dash,  color = colours[6])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n7), linestyle=:dash,  color = colours[7])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n8), linestyle=:dash,  color = colours[8])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n9), linestyle=:dash,  color = colours[9])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n10), linestyle=:dash,  color = colours[10])

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
axislegend(ax_ΔT)
#axislegend(ax_T)

fig

frames = 1:length(times)

save(end_name * "_compare.png", fig)
#save("match_Feltham_no_collisions.png", fig) 
