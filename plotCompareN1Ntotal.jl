using CairoMakie
using JLD2
using PerceptualColourMaps
using Oceananigans
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat # these constants can be called inside any function
using Statistics
include("test_functions.jl")
# change to number density * crystal volume
# turn off rise

# compare the number of n1 and the total number on two separate plots for the 6 options. Then also compare for the different choices of epsilon

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

feltham_name = "match_Feltham_no_n_max_eps_0.00001"#"match_Feltham_no_n_max_redistribute"
tmin = 0
tmax = 3

fig = Figure(size = (850, 850))

ax_ΔT = Axis(fig[1, 1:2];
            ylabel = "ΔT (ᵒC)",
            xlabel = "time (days)",
            limits = ((tmin, tmax), nothing))

ax_n1 = Axis(fig[1, 3:4];
            ylabel = "n₁",
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

#colours = cmap("Gouldian", N = length(Rs))
colours = cmap("CBTL1", N = Int(4)) 
#cvelnum = 1
#pnum = 1

for cvelnum = [1, 2, 3]
    for pnum = [1, 2]
        end_name = "match_Feltham_no_n_max_eps_00001_depth_1"
        new_label = ""

        collision_velocity_parameterisation_num = cvelnum # 1 is old cylinder, 2 is new cylinder, 3 is new spherical
        concentration_parameterisation_num = pnum # 1 is mean n, 2 is sum over nj
        crystal_size_collision_redistribution = 1 # 1 is old redistribution, 2 is new redistribution
        effective_radius = true # add in their effective radius

        if collision_velocity_parameterisation_num == 2
            new_label = "new cylinder"
            end_name = end_name * "_new_cyl"
        elseif collision_velocity_parameterisation_num == 3
            new_label = "new spherical"
            end_name = end_name * "_new_spherical"
        else
            new_label = "old cylinder"
        end
        if concentration_parameterisation_num == 2
            new_label = new_label * ", sum nj"
            end_name = end_name * "_sum_nj"
        end

        time_series = (;
            w = FieldTimeSeries("1D_fields"* end_name * ".jld2", "w"),
            T = FieldTimeSeries("1D_fields"* end_name * ".jld2", "T"),
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
        #ΔS =  time_series.S[1, 1, :, :] .- time_series.S[1, 1, :, 1]

        tmax = maximum(times)
        tmin = 0

        ntotal = find_z_mean(time_series.n₁) .+ find_z_mean(time_series.n₂) + find_z_mean(time_series.n₃) .+ find_z_mean(time_series.n4) .+ find_z_mean(time_series.n5) .+ find_z_mean(time_series.n6) .+ find_z_mean(time_series.n7) .+ find_z_mean(time_series.n8) .+ find_z_mean(time_series.n9) .+ find_z_mean(time_series.n10)

        

        if pnum == 1
            lines!(ax_ΔT, times, find_z_mean(time_series.T) .- find_z_mean(time_series.T[:, :, :, 1]), linestyle=:dash, color = colours[cvelnum], label = new_label)
            lines!(ax_T, times, find_z_mean(time_series.T) .- Tf, color = colours[cvelnum], linestyle=:dash, label = new_label)
            lines!(ax_n, times, ntotal, linestyle=:dash, color = colours[cvelnum], label = new_label)
            lines!(ax_n1, times, find_z_mean(time_series.n₁), linestyle=:dash, label = new_label, color = colours[cvelnum])
        else
            lines!(ax_n, times, ntotal, color = colours[cvelnum], label = new_label)
            lines!(ax_n1, times, find_z_mean(time_series.n₁), label = new_label, color = colours[cvelnum])
            lines!(ax_ΔT, times, find_z_mean(time_series.T) .- find_z_mean(time_series.T[:, :, :, 1]), color = colours[cvelnum], label = new_label)
            lines!(ax_T, times, find_z_mean(time_series.T) .- Tf, color = colours[cvelnum], label = new_label)
        end
        





        
    end
end

fig      
save("eps0.00001_compare_nomax_d1.png", fig)

axislegend(ax_n1, position=(:right, :bottom, ))
axislegend(ax_n, position=(:right, :bottom, ))
axislegend(ax_ΔT, position=(:right, :bottom, ))
axislegend(ax_T, position=(:right, :bottom, ))

save("eps0.00001_compare_nomax2_d1.png", fig)   