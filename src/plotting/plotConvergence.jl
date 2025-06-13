# load in the data and find the time to reach 10% of the maximum number in ...
# save the distribution and the time

using Oceananigans
using Statistics
using JLD2
using CairoMakie
using PerceptualColourMaps
using KernelDensity
include("loadData.jl")
using .loadData

include("plottingFunctions.jl")
using .plottingFunctions

colours = cmap("CBTL1", N=10) 

fig = Figure(size = (850, 600))
ax = Axis(fig[1, 1]; xlabel = "N", ylabel = "t(s)")

fig2 = Figure(size = (850, 600))
ax2 = Axis(fig2[1, 1]; xlabel = "N", ylabel = "r̅(mm)")

fig3 = Figure(size = (850, 600))
ax3 = Axis(fig3[1, 1]; xlabel = "N", ylabel = "n1")

ϵ = 1e-8
for numSizeClasses = [50, 100, 150, 200, 250]

Rs = range(0.01, 2, numSizeClasses) .* 1e-3
aspect_ratio = 50
Volume = 1

for crystal_size_collision_redistribution = [2]
for cvelnum = [3]
for pnum = [1]
    end_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
    new_label = ""

collision_velocity_parameterisation_num = cvelnum # 1 is old cylinder, 2 is new cylinder, 3 is new spherical
concentration_parameterisation_num = pnum # 1 is mean n, 2 is sum over nj
#crystal_size_collision_redistribution = 1 # 1 is old redistribution, 2 is new redistribution
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
if crystal_size_collision_redistribution == 2
    end_name = end_name * "_redistribute"
end

end_filename = "1D_fields" * end_name * ".jld2"
println(cvelnum)
println(end_name)
try
    @load end_name * "timeAndNs.jld2" timeAndNs
    meanrad2 = sum(Rs[2:end] .* timeAndNs[2:end-1])/sum(timeAndNs[2:end-1]) # mean radius from 2 upwards

    scatter!(ax3, numSizeClasses, timeAndNs[1], color = colours[pnum])
    
    if crystal_size_collision_redistribution == 2
        if cvelnum == 1
            scatter!(ax, numSizeClasses, timeAndNs[end], color = colours[3])
            scatter!(ax2, numSizeClasses, meanrad2, color = colours[3])
        elseif cvelnum == 2
            scatter!(ax, numSizeClasses, timeAndNs[end], marker = :rect, color = colours[3])
            scatter!(ax2, numSizeClasses, meanrad2, marker = :rect, color = colours[3])
        else
            scatter!(ax, numSizeClasses, timeAndNs[end], marker = :star5, color = colours[3])
            scatter!(ax2, numSizeClasses, meanrad2, marker = :star5, color = colours[3])
        end
    else
        if cvelnum == 1
            scatter!(ax, numSizeClasses, timeAndNs[end], color = colours[pnum])
            scatter!(ax2, numSizeClasses, meanrad2, color = colours[pnum])
        elseif cvelnum == 2
            scatter!(ax, numSizeClasses, timeAndNs[end], marker = :rect, color = colours[pnum])
            scatter!(ax2, numSizeClasses, meanrad2, marker = :rect, color = colours[pnum])
        else
            scatter!(ax, numSizeClasses, timeAndNs[end], marker = :star5, color = colours[pnum])
            scatter!(ax2, numSizeClasses, meanrad2, marker = :star5, color = colours[pnum])
        end
    end

catch
end

end
end
end
end

display(fig)
display(fig2)
display(fig3)
save("convergence_times.png", fig)
save("convergence_meanrad.png", fig2)
save("convergence_n1.png", fig3)

ϵ = 1e-8
for numSizeClasses = [50]#, 100, 150, 200, 250]#[1e-2, 1e-3, 1e-4, 1e-5, 1e-8]

Rs = range(0.01, 2, numSizeClasses) .* 1e-3
aspect_ratio = 50
Volume = 1

for cvelnum = [3]#[1, 2, 3]
for pnum = [1]#[1, 2]


end_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
new_label = ""

collision_velocity_parameterisation_num = cvelnum # 1 is old cylinder, 2 is new cylinder, 3 is new spherical
concentration_parameterisation_num = pnum # 1 is mean n, 2 is sum over nj
crystal_size_collision_redistribution = 2 # 1 is old redistribution, 2 is new redistribution
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
if crystal_size_collision_redistribution == 2
    end_name = end_name * "_redistribute"
end

end_filename = "1D_fields" * end_name * ".jld2"

try
    ts = load_plot_data(end_filename, numSizeClasses) #load_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_r_av.jld2", numSizeClasses)
    ts_times = ts.w.times

    n1max = 1.298846543545204e13 #maximum(find_z_mean(ts.n1))
    n1tofind = 0.99 * n1max
    tindx = argmin(abs.(find_z_mean(ts.n1) .- n1tofind))

    timeAndNs = zeros(numSizeClasses + 1)

    for n in range(1, stop=numSizeClasses, step=1)
        ndata = getproperty(ts, Symbol("n$n"))
        timeAndNs[n] = ndata[1, 1, 1, tindx]
    end

    timeAndNs[end] = ts_times[tindx]

    # Set tick positions and labels for every 20 bars
    #tick_pos = [1, 40, 80, 120, 160, 200]
    tick_pos = collect(range(0, stop = numSizeClasses, step = 20))
    tick_pos[1] = 1
    tick_labels = string.(round.(Rs[tick_pos]*1000, digits = 2))

    y_tick_pos = [1e10, 1e8, 1e6, 1e4, 1e2, 1e0] #, 1e-1]
    y_tick_labels = ["10¹⁰", "10⁸", "10⁶", "10⁴", "10²", "10⁰"] 

    fig = Figure(size = (1350, 850))
    ax = Axis(fig[1, 1]; xticks=(tick_pos, tick_labels), yticks=(y_tick_pos, y_tick_labels), yscale=log10, ylabel = "nᵢ", xlabel = "r(mm)", limits = ((0, numSizeClasses + 1), (1e-1, 2*maximum(timeAndNs))))
    bars = barplot!(ax, 1:numSizeClasses, timeAndNs[1:end-1], color = :gray, fillto =0.1, gap = 0, strokecolor = :black, strokewidth = 1)

    display(fig)



    save(end_name * "timeAndNs.png", fig)

    @save end_name * "timeAndNs.jld2" timeAndNs  # Saves variable `timeAndNs`
catch
end

end
end
end

