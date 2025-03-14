using CairoMakie
using JLD2
using PerceptualColourMaps
using Oceananigans
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat # these constants can be called inside any function
using Statistics
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

ϵ = 1e-8 #7.4 * 1e-6 #10^-3 # m²s⁻³
feltham_name = end_name = "_same_C_just_collisions_epsilon" * string(ϵ) #"_same_n_just_collisions_epsilon" * string(ϵ)


#colours = cmap("Gouldian", N = length(Rs))
colours = cmap("CBTL1", N = Int(1.2*length(Rs))) 

#new_label = "ϵ = 0.005"

feltham_time_series = (;
     w = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "w"),
     u = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "u"),
     v = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "v"),
     T = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "T"),
     S = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "S"),
     n1 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n1"),
     n2 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n2"),
     n3 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n3"),
     n4 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n4"),
     n5 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n5"),
     n6 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n6"),
     n7 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n7"),
     n8 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n8"),
     n9 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n9"),
     n10 = FieldTimeSeries("1D_fields"* feltham_name * ".jld2", "n10"),
     )

feltham_times = feltham_time_series.w.times
#feltham_times = feltham_times/60
feltham_times = feltham_times/(3600*24)

cvelnum = 1
pnum = 2
#for cvelnum = [1]#[1, 2, 3]
#for pnum = [2]#[1, 2]
end_name = "_same_C_just_collisions_epsilon" * string(ϵ)
new_label = ""

collision_velocity_parameterisation_num = cvelnum # 1 is old cylinder, 2 is new cylinder, 3 is new spherical
concentration_parameterisation_num = pnum # 1 is mean n, 2 is sum over nj
crystal_size_collision_redistribution = 1 # 1 is old redistribution, 2 is new redistribution
effective_radius = true # add in their effective radius

if collision_velocity_parameterisation_num == 2
    end_name = end_name * "_new_cyl"
    new_label = "new cylinder"
elseif collision_velocity_parameterisation_num == 3
    end_name = end_name * "_new_spherical"
        new_label = "new spherical"
end
if concentration_parameterisation_num == 2
    end_name = end_name * "_sum_nj"
    new_label = new_label * ", sum nj"
end

time_series = (;
    w = FieldTimeSeries("1D_fields"* end_name * ".jld2", "w"),
    u = FieldTimeSeries("1D_fields"* end_name * ".jld2", "u"),
    v = FieldTimeSeries("1D_fields"* end_name * ".jld2", "v"),
    T = FieldTimeSeries("1D_fields"* end_name * ".jld2", "T"),
    S = FieldTimeSeries("1D_fields"* end_name * ".jld2", "S"),
    n1 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n1"),
    n2 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n2"),
    n3 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n3"),
    n4 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n4"),
    n5 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n5"),
    n6 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n6"),
    n7 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n7"),
    n8 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n8"),
    n9 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n9"),
    n10 = FieldTimeSeries("1D_fields"* end_name * ".jld2", "n10"),
    )

times = time_series.w.times
#times = times/60
times = times/(3600*24)


tmax =  maximum(times)
tmin = 0

fig = Figure(size = (850, 850))

ax_no = Axis(fig[2, 1:2];
            ylabel = "n, mean density",
            xlabel = "time (days)",
            yscale = log10,
            #limits = ((0.01, 0.02), (1e4, 1e6)))
            limits = ((tmin, tmax), (0.1, 1e13)))

ax_C = Axis(fig[1, 3:4];
            ylabel = "C",
            xlabel = "time (days)",
            yscale = log10,
            limits = ((tmin, tmax), (1e-15, 2e-3)))

ax_Co = Axis(fig[1, 1:2],
            ylabel = "C, mean density",
            xlabel = "time (days)",
            yscale = log10,
            limits = ((tmin, tmax), (1e-15, 2e-3)))

ax_n = Axis(fig[2, 3:4];
            ylabel = "n",
            xlabel = "time (days)",
            yscale = log10,
            #limits = ((tmin, tmax), (1e4, 1e6)))
            limits = ((tmin, tmax), (0.1, 1e13)))


lines!(ax_Co, feltham_times, find_C_const(1)*find_z_mean(feltham_time_series.n1), linestyle=:solid,  color = colours[1])
lines!(ax_Co, feltham_times, find_C_const(2)*find_z_mean(feltham_time_series.n2), linestyle=:solid,  color = colours[2])
lines!(ax_Co, feltham_times, find_C_const(3)*find_z_mean(feltham_time_series.n3), linestyle=:solid,  color = colours[3])
lines!(ax_Co, feltham_times, find_C_const(4)*find_z_mean(feltham_time_series.n4), linestyle=:solid,  color = colours[4])
lines!(ax_Co, feltham_times, find_C_const(5)*find_z_mean(feltham_time_series.n5), linestyle=:solid,  color = colours[5])
lines!(ax_Co, feltham_times, find_C_const(6)*find_z_mean(feltham_time_series.n6), linestyle=:solid,  color = colours[6])
lines!(ax_Co, feltham_times, find_C_const(7)*find_z_mean(feltham_time_series.n7), linestyle=:solid,  color = colours[7])
lines!(ax_Co, feltham_times, find_C_const(8)*find_z_mean(feltham_time_series.n8), linestyle=:solid,  color = colours[8])
lines!(ax_Co, feltham_times, find_C_const(9)*find_z_mean(feltham_time_series.n9), linestyle=:solid,  color = colours[9])
lines!(ax_Co, feltham_times, find_C_const(10)*find_z_mean(feltham_time_series.n10), linestyle=:solid, color = colours[10])

lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n1), linestyle=:solid, color = colours[1])
lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n2), linestyle=:solid,  color = colours[2])
lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n3), linestyle=:solid,  color = colours[3])
lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n4), linestyle=:solid,  color = colours[4])
lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n5), linestyle=:solid,  color = colours[5])
lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n6), linestyle=:solid,  color = colours[6])
lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n7), linestyle=:solid,  color = colours[7])
lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n8), linestyle=:solid,  color = colours[8])
lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n9), linestyle=:solid,  color = colours[9])
lines!(ax_no, feltham_times, find_z_mean(feltham_time_series.n10), linestyle=:solid,  color = colours[10])

lines!(ax_C, times, find_C_const(1)*find_z_mean(time_series.n1), label = "C₁", color = colours[1])
lines!(ax_C, times, find_C_const(2)*find_z_mean(time_series.n2), label = "C₂", color = colours[2])
lines!(ax_C, times, find_C_const(3)*find_z_mean(time_series.n3), label = "C₃", color = colours[3])
lines!(ax_C, times, find_C_const(4)*find_z_mean(time_series.n4), label = "C₄", color = colours[4])
lines!(ax_C, times, find_C_const(5)*find_z_mean(time_series.n5), label = "C₅", color = colours[5])
lines!(ax_C, times, find_C_const(6)*find_z_mean(time_series.n6), label = "C₆", color = colours[6])
lines!(ax_C, times, find_C_const(7)*find_z_mean(time_series.n7), label = "C₇", color = colours[7])
lines!(ax_C, times, find_C_const(8)*find_z_mean(time_series.n8), label = "C₈", color = colours[8])
lines!(ax_C, times, find_C_const(9)*find_z_mean(time_series.n9), label = "C₉", color = colours[9])
lines!(ax_C, times, find_C_const(10)*find_z_mean(time_series.n10), label = "C₁₀", color = colours[10])

lines!(ax_n, times, find_z_mean(time_series.n1), label = "n₁", color = colours[1])
lines!(ax_n, times, find_z_mean(time_series.n2), label = "n₂", color = colours[2])
lines!(ax_n, times, find_z_mean(time_series.n3), label = "n₃", color = colours[3])
lines!(ax_n, times, find_z_mean(time_series.n4), label = "n₄", color = colours[4])
lines!(ax_n, times, find_z_mean(time_series.n5), label = "n₅", color = colours[5])
lines!(ax_n, times, find_z_mean(time_series.n6), label = "n₆", color = colours[6])
lines!(ax_n, times, find_z_mean(time_series.n7), label = "n₇", color = colours[7])
lines!(ax_n, times, find_z_mean(time_series.n8), label = "n₈", color = colours[8])
lines!(ax_n, times, find_z_mean(time_series.n9), label = "n₉", color = colours[9])
lines!(ax_n, times, find_z_mean(time_series.n10), label = "n₁₀", color = colours[10])


axislegend(ax_C)
axislegend(ax_n)

fig

frames = 1:length(times)

#save(end_name * "_4_compare_zoom_max.png", fig)
save(end_name * "_compare_nomax.png", fig)
#save("match_Feltham_no_collisions.png", fig) 

#end
#end

# plot the figures on top of each other


# also plot the mean number and mean concentration
# plot mean radius against time (one line including smallest size class, other line excluding smallest size class)

# Efficient computation of weighted sum
mean_radius_f = sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) .* Rs[i] for i in 1:10)./sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) for i in 1:10)
mean_radius_exclude_small_f = sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) .* Rs[i] for i in 2:10)./sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) for i in 2:10)
num_g1_f = sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) for i in 2:10)

mean_radius = sum(find_z_mean(getfield(time_series, Symbol("n$i"))) .* Rs[i] for i in 1:10)./sum(find_z_mean(getfield(time_series, Symbol("n$i"))) for i in 1:10)
mean_radius_exclude_small = sum(find_z_mean(getfield(time_series, Symbol("n$i"))) .* Rs[i] for i in 2:10)./sum(find_z_mean(getfield(time_series, Symbol("n$i"))) for i in 2:10)
num_g1 = sum(find_z_mean(getfield(time_series, Symbol("n$i"))) for i in 2:10)

fig2 = Figure(size = (850, 450))

ax_ro = Axis(fig2[1, 1];
            ylabel ="r̅",
            xlabel = "time (minutes)",
            yscale = log10,
            limits = ((tmin, tmax), (9e-6, 1e-2)))

ax_r = Axis(fig2[1, 2];
            ylabel = "r̅",
            xlabel = "time (minutes)",
            yscale = log10,
            limits = ((tmin, tmax), (9e-6, 1e-2)))


lines!(ax_ro, feltham_times, mean_radius_f, linestyle=:solid,  color = colours[1], label = "all i")
lines!(ax_ro, feltham_times[num_g1_f .> 0], mean_radius_exclude_small_f[num_g1_f .> 0], linestyle=:solid,  color = colours[9], label = "i > 1")

lines!(ax_r, times, mean_radius, linestyle=:solid,  color = colours[1], label = "all i")
lines!(ax_r, times[num_g1 .> 0], mean_radius_exclude_small[num_g1 .> 0], linestyle=:solid,  color = colours[9], label = "i > 1")

axislegend(ax_ro)

fig2

frames = 1:length(times)

save(end_name * "_compare_rs_zoom.png", fig2)

conc = sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) .* 2*π*Rs[i]^3/50 for i in 1:10)

conc = sum(find_z_mean(getfield(time_series, Symbol("n$i"))) .* 2*π*Rs[i]^3/50 for i in 1:10)

#################################################

fig = Figure(size = (1050, 450))

ax_r = Axis(fig[1, 5:6];
            ylabel ="r̅",
            xlabel = "time (minutes)",
            yscale = log10,
            limits = ((tmin, tmax), (9e-6, 1e-2)))

ax_C = Axis(fig[1, 3:4];
            ylabel = "C",
            xlabel = "time (minutes)",
            yscale = log10,
            limits = ((tmin, tmax), (1e-15, 2e-3)))

ax_n = Axis(fig[1, 1:2];
            ylabel = "n",
            xlabel = "time (minutes)",
            yscale = log10,
            #limits = ((tmin, tmax), (1e4, 1e6)))
            limits = ((tmin, tmax), (0.1, 5e13)))


lines!(ax_C, feltham_times, find_C_const(1)*find_z_mean(feltham_time_series.n1), linestyle=:dash,  color = colours[1])
lines!(ax_C, feltham_times, find_C_const(2)*find_z_mean(feltham_time_series.n2), linestyle=:dash,  color = colours[2])
lines!(ax_C, feltham_times, find_C_const(3)*find_z_mean(feltham_time_series.n3), linestyle=:dash,  color = colours[3])
lines!(ax_C, feltham_times, find_C_const(4)*find_z_mean(feltham_time_series.n4), linestyle=:dash,  color = colours[4])
lines!(ax_C, feltham_times, find_C_const(5)*find_z_mean(feltham_time_series.n5), linestyle=:dash,  color = colours[5])
lines!(ax_C, feltham_times, find_C_const(6)*find_z_mean(feltham_time_series.n6), linestyle=:dash,  color = colours[6])
lines!(ax_C, feltham_times, find_C_const(7)*find_z_mean(feltham_time_series.n7), linestyle=:dash,  color = colours[7])
lines!(ax_C, feltham_times, find_C_const(8)*find_z_mean(feltham_time_series.n8), linestyle=:dash,  color = colours[8])
lines!(ax_C, feltham_times, find_C_const(9)*find_z_mean(feltham_time_series.n9), linestyle=:dash,  color = colours[9])
lines!(ax_C, feltham_times, find_C_const(10)*find_z_mean(feltham_time_series.n10), linestyle=:dash, color = colours[10])

lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n1), linestyle=:dash, color = colours[1])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n2), linestyle=:dash,  color = colours[2])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n3), linestyle=:dash,  color = colours[3])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n4), linestyle=:dash,  color = colours[4])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n5), linestyle=:dash,  color = colours[5])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n6), linestyle=:dash,  color = colours[6])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n7), linestyle=:dash,  color = colours[7])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n8), linestyle=:dash,  color = colours[8])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n9), linestyle=:dash,  color = colours[9])
lines!(ax_n, feltham_times, find_z_mean(feltham_time_series.n10), linestyle=:dash,  color = colours[10])

lines!(ax_C, times, find_C_const(1)*find_z_mean(time_series.n1), label = "C₁", color = colours[1])
lines!(ax_C, times, find_C_const(2)*find_z_mean(time_series.n2), label = "C₂", color = colours[2])
lines!(ax_C, times, find_C_const(3)*find_z_mean(time_series.n3), label = "C₃", color = colours[3])
lines!(ax_C, times, find_C_const(4)*find_z_mean(time_series.n4), label = "C₄", color = colours[4])
lines!(ax_C, times, find_C_const(5)*find_z_mean(time_series.n5), label = "C₅", color = colours[5])
lines!(ax_C, times, find_C_const(6)*find_z_mean(time_series.n6), label = "C₆", color = colours[6])
lines!(ax_C, times, find_C_const(7)*find_z_mean(time_series.n7), label = "C₇", color = colours[7])
lines!(ax_C, times, find_C_const(8)*find_z_mean(time_series.n8), label = "C₈", color = colours[8])
lines!(ax_C, times, find_C_const(9)*find_z_mean(time_series.n9), label = "C₉", color = colours[9])
lines!(ax_C, times, find_C_const(10)*find_z_mean(time_series.n10), label = "C₁₀", color = colours[10])

lines!(ax_n, times, find_z_mean(time_series.n1), label = "n₁", color = colours[1])
lines!(ax_n, times, find_z_mean(time_series.n2), label = "n₂", color = colours[2])
lines!(ax_n, times, find_z_mean(time_series.n3), label = "n₃", color = colours[3])
lines!(ax_n, times, find_z_mean(time_series.n4), label = "n₄", color = colours[4])
lines!(ax_n, times, find_z_mean(time_series.n5), label = "n₅", color = colours[5])
lines!(ax_n, times, find_z_mean(time_series.n6), label = "n₆", color = colours[6])
lines!(ax_n, times, find_z_mean(time_series.n7), label = "n₇", color = colours[7])
lines!(ax_n, times, find_z_mean(time_series.n8), label = "n₈", color = colours[8])
lines!(ax_n, times, find_z_mean(time_series.n9), label = "n₉", color = colours[9])
lines!(ax_n, times, find_z_mean(time_series.n10), label = "n₁₀", color = colours[10])

lines!(ax_r, feltham_times, mean_radius_f, linestyle=:dash,  color = colours[1])
lines!(ax_r, feltham_times[num_g1_f .> 0], mean_radius_exclude_small_f[num_g1_f .> 0], linestyle=:dash,  color = colours[9])

lines!(ax_r, times, mean_radius, linestyle=:solid,  color = colours[1], label = "all i")
lines!(ax_r, times[num_g1 .> 0], mean_radius_exclude_small[num_g1 .> 0], linestyle=:solid,  color = colours[9], label = "i > 1")

axislegend(ax_r)


axislegend(ax_C)
axislegend(ax_n)

fig

frames = 1:length(times)

#save(end_name * "_4_compare_zoom_max.png", fig)
save(end_name * "_compare_nomax_ontop.png", fig)
