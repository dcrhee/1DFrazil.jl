using CairoMakie
using JLD2
using Oceananigans
using PerceptualColourMaps

# change to number density * crystal volume
# turn off rise

Uₐ =  5 # m/s specify the wind strength
Fetch = 1500 # m wind fetch
Tₐ = -20 # atmosphere temperature

end_name = "_Ta_" * string(Tₐ) * "_X_" * string(Fetch) * "_Ua_" * string(Uₐ)
end_name = "no_vel"

end_name = "match_Feltham_no_n_max_eps_0.00001_depth_10_1D"

numz = 128

Rs = [0.01, 0.05, 0.15, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 2] .* 1e-3
colours = cmap("CBTL1", N = Int(1.2*length(Rs))) 

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
    n4 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n4"),
    n5 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n5"),
    n6 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n6"),
    n7 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n7"),
    n8 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n8"),
    n9 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n9"),
    n10 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n10"),
)

times = time_series.w.times

#ΔT = time_series.T[1, 1, :, :] .- time_series.T[1, 1, :, 1]
#ΔS =  time_series.S[1, 1, :, :] .- time_series.S[1, 1, :, 1]

xw, yw, zw = nodes(time_series.w)
xu, yu, zu = nodes(time_series.u)
xv, yv, zv = nodes(time_series.v)
xT, yT, zT = nodes(time_series.T)
xS, yS, zS = nodes(time_series.S)
xC, yC, zC = nodes(time_series.n₁)

n = Observable(1)

fig = Figure(size = (850, 850))

ax_ΔT = Axis(fig[1, 1:2];
            xlabel = "ΔT (ᵒC)",
            ylabel = "z (m)",)
            #limits = ((minimum(ΔT), maximum(ΔT)), nothing))

#ax_ΔS = Axis(fig[1, 3:4];
#            xlabel = "ΔS (ppt)",
#            ylabel = "z (m)",
#            limits = ((minimum(ΔS), maximum(ΔS)), nothing))

ax_w = Axis(fig[1, 5:6];
            xlabel = "w (m s⁻¹)",
            ylabel = "z (m)",)
            #limits = ((minimum(time_series.w), maximum(time_series.w)), nothing))

ax_T = Axis(fig[2, 1:2];
              xlabel = "T (ᵒC)",
              ylabel = "z (m)",)
#              limits = ((minimum(time_series.T), maximum(time_series.T)), nothing))

ax_S = Axis(fig[2, 3:4];
              xlabel = "S (psu)",
              ylabel = "z (m)",)
              #limits = ((minimum(time_series.S), maximum(time_series.S)), nothing))

ax_C = Axis(fig[2, 5:6];
              xlabel = "C",
              ylabel = "z (m)",)
#              limits = ((minimum(time_series.n₁), maximum(time_series.n₁)), nothing))


uₙ = @lift time_series.u[$n][1, 1, 1:numz]
vₙ = @lift time_series.v[$n][1, 1, 1:numz]
wₙ = @lift time_series.w[$n][1, 1, 1:numz]
Tₙ = @lift time_series.T[$n][1, 1, 1:numz]
Sₙ = @lift time_series.S[$n][1, 1, 1:numz]
n₁ₙ = @lift time_series.n₁[$n][1, 1, 1:numz]
n₂ₙ = @lift time_series.n₂[$n][1, 1, 1:numz]
n₃ₙ = @lift time_series.n₃[$n][1, 1, 1:numz]
n4ₙ = @lift time_series.n4[$n][1, 1, 1:numz]
n5ₙ = @lift time_series.n5[$n][1, 1, 1:numz]
n6ₙ = @lift time_series.n6[$n][1, 1, 1:numz]
n7ₙ = @lift time_series.n7[$n][1, 1, 1:numz]
n8ₙ = @lift time_series.n8[$n][1, 1, 1:numz]
n9ₙ = @lift time_series.n9[$n][1, 1, 1:numz]
n10ₙ = @lift time_series.n10[$n][1, 1, 1:numz]
ΔTₙ = @lift time_series.T[$n][1, 1, 1:numz] - time_series.T[1][1, 1, 1:numz]
ΔSₙ = @lift time_series.S[$n][1, 1, 1:numz] - time_series.S[1][1, 1, 1:numz]

lines!(ax_ΔT, ΔTₙ, zT)
lines!(ax_ΔS, ΔSₙ, zS)
#lines!(ax_w, wₙ, zw)
lines!(ax_T, Tₙ, zT)
lines!(ax_S, Sₙ, zS)
lines!(ax_C, n₁ₙ, zC, label = "n₁", color = colours[1])
lines!(ax_C, n₂ₙ, zC, label = "n₂", color = colours[2])
lines!(ax_C, n₃ₙ, zC, label = "n₃", color = colours[3])
lines!(ax_C, n4ₙ, zC, label = "n4", color = colours[4])
lines!(ax_C, n5ₙ, zC, label = "n5", color = colours[5])
lines!(ax_C, n6ₙ, zC, label = "n6", color = colours[6])
lines!(ax_C, n7ₙ, zC, label = "n7", color = colours[7])
lines!(ax_C, n8ₙ, zC, label = "n8", color = colours[8])
lines!(ax_C, n9ₙ, zC, label = "n9", color = colours[9])
lines!(ax_C, n10ₙ, zC, label = "n10", color = colours[10])


axislegend(ax_C)

fig

frames = 1:length(times)

record(fig, "rise_noflux10.mp4", frames, framerate=8) do i
    n[] = i
end