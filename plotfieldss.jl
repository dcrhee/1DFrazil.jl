using CairoMakie
using JLD2

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
     n₃ = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n₃"))

times = time_series.w.times

xw, yw, zw = nodes(time_series.w)
xu, yu, zu = nodes(time_series.u)
xv, yv, zv = nodes(time_series.v)
xT, yT, zT = nodes(time_series.T)
xS, yS, zS = nodes(time_series.S)
xC, yC, zC = nodes(time_series.n₁)

n = Observable(1)

fig = Figure(size = (850, 850))

ax_u = Axis(fig[1, 1:2];
            xlabel = "u (m s⁻¹)",
            ylabel = "z (m)",
            limits = ((minimum(time_series.u), maximum(time_series.u)), nothing))

ax_v = Axis(fig[1, 3:4];
            xlabel = "v (m s⁻¹)",
            ylabel = "z (m)",
            limits = ((minimum(time_series.v), maximum(time_series.v)), nothing))

ax_w = Axis(fig[1, 5:6];
            xlabel = "w (m s⁻¹)",
            ylabel = "z (m)",
            limits = ((minimum(time_series.w), maximum(time_series.w)), nothing))

ax_T = Axis(fig[2, 1:2];
              xlabel = "T (ᵒC)",
              ylabel = "z (m)",
              limits = ((minimum(time_series.T), maximum(time_series.T)), nothing))

ax_S = Axis(fig[2, 3:4];
              xlabel = "S (psu)",
              ylabel = "z (m)",)
              #limits = ((minimum(time_series.S), maximum(time_series.S)), nothing))

ax_C = Axis(fig[2, 5:6];
              xlabel = "C",
              ylabel = "z (m)",
              limits = ((minimum(time_series.n₂), maximum(time_series.n₂)), nothing))

uₙ = @lift time_series.u[$n][1, 1, :]
vₙ = @lift time_series.v[$n][1, 1, :]
wₙ = @lift time_series.w[$n][1, 1, :]
Tₙ = @lift time_series.T[$n][1, 1, :]
Sₙ = @lift time_series.S[$n][1, 1, :]
n₁ₙ = @lift time_series.n₁[$n][1, 1, :]
n₂ₙ = @lift time_series.n₂[$n][1, 1, :]
n₃ₙ = @lift time_series.n₃[$n][1, 1, :]

lines!(ax_u, uₙ, zu)
lines!(ax_v, vₙ, zv)
lines!(ax_w, wₙ, zw)
lines!(ax_T, Tₙ, zT)
lines!(ax_S, Sₙ, zS)
lines!(ax_C, n₁ₙ, zC, label = "n₁")
lines!(ax_C, n₂ₙ, zC, label = "n₂")
lines!(ax_C, n₃ₙ, zC, label = "n₃")

axislegend(ax_C)

fig

frames = 1:length(times)

record(fig, "test.mp4", frames, framerate=8) do i
    n[] = i
end