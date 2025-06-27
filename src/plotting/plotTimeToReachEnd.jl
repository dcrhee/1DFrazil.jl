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

colours = cmap("CBTL1", N=8) 

fig = Figure(size = (850, 600))
ax = Axis(fig[1, 1]; xscale=log10, xlabel = "ϵ (m²s⁻³)", ylabel = "t (s)")

fig2 = Figure(size = (850, 600))
ax2 = Axis(fig2[1, 1]; xscale=log10, xlabel = "ϵ (m²s⁻³)", ylabel = "r̅ (mm)")

fig3 = Figure(size = (850, 600))
ax3 = Axis(fig3[1, 1]; xscale=log10, xlabel = "ϵ (m²s⁻³)", ylabel = "n1")

fig = Figure(size = (850, 600))
ax = Axis(fig[1, 1]; xscale=log10, xlabel = "ϵ (m²s⁻³)", ylabel = "t (s)")
ax2 = Axis(fig[1, 2]; xscale=log10, xlabel = "ϵ (m²s⁻³)", ylabel = "r̅ (mm)")


msize = 10

p1c1_times = []
p1c1_epsilon = []
p1c1_rad = []

p2c3_times = []
p2c3_epsilon = []
p2c3_rad = []

for ϵ = [1e-2, 1e-3, 1e-4, 1e-5, 1e-8]

numSizeClasses = 200
#Rs = range(0.01, 2, numSizeClasses) .* 1e-3
Rs = exp.(range(start=log(0.01*0.001), stop=log(0.002), length=numSizeClasses))
aspect_ratio = 50
Volume = 1

for crystal_size_collision_redistribution = [1] #, 2]
for cvelnum = [1, 2, 3]
for pnum = [1, 2]
    if pnum == 1
        if crystal_size_collision_redistribution == 1
            cindx = 1
        else
            cindx = 4
        end
    else
        if crystal_size_collision_redistribution == 1
            cindx = 6
        else
            cindx = 5
        end
    end
    #end_name = "sameConc_fast_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
    end_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
    #end_name = "log_space_sameConc_fast_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
    new_label = ""
    label_str = ""

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
else
    new_label = "old cylinder"
end
if concentration_parameterisation_num == 2
    end_name = end_name * "_sum_nj"
    new_label = new_label * ", differential crystal size"
else
    new_label = new_label * ", point-like"
end
if crystal_size_collision_redistribution == 2
    end_name = end_name * "_redistribute"
end

end_filename = "1D_fields" * end_name * ".jld2"
println(cvelnum)
println(end_name)
try
    @load end_name * "timeAndNs.jld2" timeAndNs
    meanrad2 = sum(Rs[2:end] .* timeAndNs[2:end-1])/sum(timeAndNs[2:end-1])*1000 # mean radius from 2 upwards

    if crystal_size_collision_redistribution == 1
        if pnum == 1
            if cvelnum == 1
                push!(p1c1_times, timeAndNs[end])
                push!(p1c1_epsilon, ϵ)
                push!(p1c1_rad, meanrad2)
            end
        else
            if cvelnum == 3
                push!(p2c3_times, timeAndNs[end])
                push!(p2c3_epsilon, ϵ)
                push!(p2c3_rad, meanrad2)
            end

        end
    end
    
    scatter!(ax3, ϵ, timeAndNs[1], color = colours[cindx])
    if ϵ == 0.01
        if crystal_size_collision_redistribution == 2
            if cvelnum == 1
                scatter!(ax, ϵ, timeAndNs[end], marker = :diamond, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax4, ϵ, meanrad2, marker = :diamond, color = colours[cindx], markersize = msize, label = new_label)
            elseif cvelnum == 2
                scatter!(ax, ϵ, timeAndNs[end], marker = :star8, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax4, ϵ, meanrad2, marker = :star8, color = colours[cindx], markersize = msize, label = new_label)
            else
                scatter!(ax, ϵ, timeAndNs[end], marker = :xcross, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax4, ϵ, meanrad2, marker = :xcross, color = colours[cindx], markersize = msize, label = new_label)
            end
        else
            if cvelnum == 1
                scatter!(ax, ϵ, timeAndNs[end], color = colours[cindx], marker = :rect, markersize = msize, label = new_label)
                scatter!(ax2, ϵ, meanrad2, color = colours[cindx], marker = :rect, markersize = msize, label = new_label)
            elseif cvelnum == 2
                scatter!(ax, ϵ, timeAndNs[end], marker = :circle, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax2, ϵ, meanrad2, marker = :circle, color = colours[cindx], markersize = msize, label = new_label)
            else
                scatter!(ax, ϵ, timeAndNs[end], marker = :cross, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax2, ϵ, meanrad2, marker = :cross, color = colours[cindx], markersize = msize, label = new_label)
            end
        end
    else
        if crystal_size_collision_redistribution == 2
            if cvelnum == 1
                scatter!(ax, ϵ, timeAndNs[end], marker = :diamond, color = colours[cindx], markersize = msize)
                scatter!(ax4, ϵ, meanrad2, marker = :diamond, color = colours[cindx], markersize = msize)
            elseif cvelnum == 2
                scatter!(ax, ϵ, timeAndNs[end], marker = :star8, color = colours[cindx], markersize = msize)
                scatter!(ax4, ϵ, meanrad2, marker = :star8, color = colours[cindx], markersize = msize)
            else
                scatter!(ax, ϵ, timeAndNs[end], marker = :xcross, color = colours[cindx], markersize = msize)
                scatter!(ax4, ϵ, meanrad2, marker = :xcross, color = colours[cindx], markersize = msize)
            end
        else
            if cvelnum == 1
                scatter!(ax, ϵ, timeAndNs[end], color = colours[cindx], marker = :rect, markersize = msize)
                scatter!(ax2, ϵ, meanrad2, color = colours[cindx], marker = :rect, markersize = msize)
            elseif cvelnum == 2
                scatter!(ax, ϵ, timeAndNs[end], marker = :circle, color = colours[cindx], markersize = msize)
                scatter!(ax2, ϵ, meanrad2, marker = :circle, color = colours[cindx], markersize = msize)
            else
                scatter!(ax, ϵ, timeAndNs[end], marker = :cross, color = colours[cindx], markersize = msize)
                scatter!(ax2, ϵ, meanrad2, marker = :cross, color = colours[cindx], markersize = msize)
            end
        end
    end

catch
end

end
end
end
end
lines!(ax, p1c1_epsilon, p1c1_times, color = colours[1])
lines!(ax, p2c3_epsilon, p2c3_times, color = colours[6])

lines!(ax2, p1c1_epsilon, p1c1_rad, color = colours[1])
lines!(ax2, p2c3_epsilon, p2c3_rad, color = colours[6])

axislegend(ax, position=(:left, :bottom), framevisible = false)
#axislegend(ax2, position=(:left, :bottom), framevisible = false)

display(fig)
#display(fig2)
#display(fig3)
#save("times.pdf", fig)
#save("meanrad.pdf", fig2)
#save("n1.pdf", fig3)

#save("times_meanrad.pdf", fig)
save("samen_times_meanrad_resitribution.pdf", fig)
#save("samen_times_meanrad_resitribution_both.pdf", fig)
#save("sameC_times_meanrad_resitribution_both.pdf", fig)
#save("logC_times_meanrad_resitribution_both.pdf", fig)

fig = Figure(size = (850, 600))
ax = Axis(fig[1:2, 1]; xscale=log10, xlabel = "ϵ (m²s⁻³)", ylabel = "t (s)")
ax2 = Axis(fig[1, 2]; xscale=log10, xlabel = "ϵ (m²s⁻³)", ylabel = "r̅ (mm)")
ax4 = Axis(fig[2, 2]; xscale=log10, xlabel = "ϵ (m²s⁻³)", ylabel = "r̅ (mm)")


msize = 10

p1c1_times = []
p1c1_epsilon = []
p1c1_rad = []

p2c3_times = []
p2c3_epsilon = []
p2c3_rad = []

for ϵ = [1e-2, 1e-3, 1e-4, 1e-5, 1e-8]

numSizeClasses = 200
#Rs = range(0.01, 2, numSizeClasses) .* 1e-3
Rs = exp.(range(start=log(0.01*0.001), stop=log(0.002), length=numSizeClasses))
aspect_ratio = 50
Volume = 1

for crystal_size_collision_redistribution = [1, 2]
for cvelnum = [1, 2, 3]
for pnum = [1, 2]
    if pnum == 1
        if crystal_size_collision_redistribution == 1
            cindx = 1
        else
            cindx = 4
        end
    else
        if crystal_size_collision_redistribution == 1
            cindx = 6
        else
            cindx = 5
        end
    end
    #end_name = "sameConc_fast_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
    #end_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
    end_name = "log_space_sameConc_fast_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
    new_label = ""
    label_str = ""

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
else
    new_label = "old cylinder"
end
if concentration_parameterisation_num == 2
    end_name = end_name * "_sum_nj"
    new_label = new_label * ", differential crystal size"
else
    new_label = new_label * ", point-like"
end
if crystal_size_collision_redistribution == 2
    end_name = end_name * "_redistribute"
end

end_filename = "1D_fields" * end_name * ".jld2"
println(cvelnum)
println(end_name)
try
    @load end_name * "timeAndNs.jld2" timeAndNs
    meanrad2 = sum(Rs[2:end] .* timeAndNs[2:end-1])/sum(timeAndNs[2:end-1])*1000 # mean radius from 2 upwards

    if crystal_size_collision_redistribution == 1
        if pnum == 1
            if cvelnum == 1
                push!(p1c1_times, timeAndNs[end])
                push!(p1c1_epsilon, ϵ)
                push!(p1c1_rad, meanrad2)
            end
        else
            if cvelnum == 3
                push!(p2c3_times, timeAndNs[end])
                push!(p2c3_epsilon, ϵ)
                push!(p2c3_rad, meanrad2)
            end

        end
    end
    
    scatter!(ax3, ϵ, timeAndNs[1], color = colours[cindx])
    if ϵ == 0.01
        if crystal_size_collision_redistribution == 2
            if cvelnum == 1
                scatter!(ax, ϵ, timeAndNs[end], marker = :diamond, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax4, ϵ, meanrad2, marker = :diamond, color = colours[cindx], markersize = msize, label = new_label)
            elseif cvelnum == 2
                scatter!(ax, ϵ, timeAndNs[end], marker = :star8, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax4, ϵ, meanrad2, marker = :star8, color = colours[cindx], markersize = msize, label = new_label)
            else
                scatter!(ax, ϵ, timeAndNs[end], marker = :xcross, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax4, ϵ, meanrad2, marker = :xcross, color = colours[cindx], markersize = msize, label = new_label)
            end
        else
            if cvelnum == 1
                scatter!(ax, ϵ, timeAndNs[end], color = colours[cindx], marker = :rect, markersize = msize, label = new_label)
                scatter!(ax2, ϵ, meanrad2, color = colours[cindx], marker = :rect, markersize = msize, label = new_label)
            elseif cvelnum == 2
                scatter!(ax, ϵ, timeAndNs[end], marker = :circle, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax2, ϵ, meanrad2, marker = :circle, color = colours[cindx], markersize = msize, label = new_label)
            else
                scatter!(ax, ϵ, timeAndNs[end], marker = :cross, color = colours[cindx], markersize = msize, label = new_label)
                scatter!(ax2, ϵ, meanrad2, marker = :cross, color = colours[cindx], markersize = msize, label = new_label)
            end
        end
    else
        if crystal_size_collision_redistribution == 2
            if cvelnum == 1
                scatter!(ax, ϵ, timeAndNs[end], marker = :diamond, color = colours[cindx], markersize = msize)
                scatter!(ax4, ϵ, meanrad2, marker = :diamond, color = colours[cindx], markersize = msize)
            elseif cvelnum == 2
                scatter!(ax, ϵ, timeAndNs[end], marker = :star8, color = colours[cindx], markersize = msize)
                scatter!(ax4, ϵ, meanrad2, marker = :star8, color = colours[cindx], markersize = msize)
            else
                scatter!(ax, ϵ, timeAndNs[end], marker = :xcross, color = colours[cindx], markersize = msize)
                scatter!(ax4, ϵ, meanrad2, marker = :xcross, color = colours[cindx], markersize = msize)
            end
        else
            if cvelnum == 1
                scatter!(ax, ϵ, timeAndNs[end], color = colours[cindx], marker = :rect, markersize = msize)
                scatter!(ax2, ϵ, meanrad2, color = colours[cindx], marker = :rect, markersize = msize)
            elseif cvelnum == 2
                scatter!(ax, ϵ, timeAndNs[end], marker = :circle, color = colours[cindx], markersize = msize)
                scatter!(ax2, ϵ, meanrad2, marker = :circle, color = colours[cindx], markersize = msize)
            else
                scatter!(ax, ϵ, timeAndNs[end], marker = :cross, color = colours[cindx], markersize = msize)
                scatter!(ax2, ϵ, meanrad2, marker = :cross, color = colours[cindx], markersize = msize)
            end
        end
    end

catch
end

end
end
end
end
lines!(ax, p1c1_epsilon, p1c1_times, color = colours[1])
lines!(ax, p2c3_epsilon, p2c3_times, color = colours[6])

lines!(ax2, p1c1_epsilon, p1c1_rad, color = colours[1])
lines!(ax2, p2c3_epsilon, p2c3_rad, color = colours[6])

axislegend(ax, position=(:left, :bottom), framevisible = false)
#axislegend(ax2, position=(:left, :bottom), framevisible = false)

display(fig)
#display(fig2)
#display(fig3)
#save("times.pdf", fig)
#save("meanrad.pdf", fig2)
#save("n1.pdf", fig3)

#save("times_meanrad.pdf", fig)
#save("samen_times_meanrad_resitribution.pdf", fig)
#save("samen_times_meanrad_resitribution_both.pdf", fig)
#save("sameC_times_meanrad_resitribution_both.pdf", fig)
save("logC_times_meanrad_resitribution_both.pdf", fig)


for ϵ = [1e-7] #, 1e-3, 1e-4, 1e-5, 1e-8]

numSizeClasses = 200
Rs = range(0.01, 2, numSizeClasses) .* 1e-3
#Rs = exp.(range(start=log(0.01*0.001), stop=log(0.002), length=numSizeClasses))
aspect_ratio = 50
Volume = 1

for crystal_size_collision_redistribution = [2]#[1, 2]
for cvelnum = [1]#[1, 2, 3]
for pnum = [2]#[1, 2]

#end_name = "sameConc_fast_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
end_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
#end_name = "fast_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
#end_name = "log_space_sameConc_fast_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
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
try
    @load end_name * "timeAndNs.jld2" timeAndNs
    #@load "sameConc_" * end_name * "timeAndNs.jld2" timeAndNs
catch
try
    println("end_name", end_name)
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
    tick_pos = [1, 40, 80, 120, 160, 200]
    tick_pos = collect(range(0, stop = numSizeClasses, step = 20))
    tick_pos[1] = 1
    tick_labels = string.(round.(Rs[tick_pos]*1000, digits = 2))

    y_tick_pos = [1e10, 1e8, 1e6, 1e4, 1e2, 1e0] #, 1e-1]
    y_tick_labels = ["10¹⁰", "10⁸", "10⁶", "10⁴", "10²", "10⁰"] 

    fig = Figure(size = (1350, 850))
    ax = Axis(fig[1, 1]; xticks=(tick_pos, tick_labels), yticks=(y_tick_pos, y_tick_labels), yscale=log10, ylabel = "nᵢ", xlabel = "r(mm)", limits = ((0, numSizeClasses + 1), (1e-1, 2*maximum(timeAndNs))))
    bars = barplot!(ax, 1:numSizeClasses, timeAndNs[1:end-1], color = :gray, fillto =0.1, gap = 0, strokecolor = :black, strokewidth = 1)

    display(fig)


    #end_name = "sameConc_" * end_name
    save(end_name * "timeAndNs.png", fig)

    @save end_name * "timeAndNs.jld2" timeAndNs  # Saves variable `timeAndNs`
catch
end
end

end
end
end
end

