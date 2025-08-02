using CairoMakie
using JLD2
using PerceptualColourMaps
using Oceananigans
include("Constants.jl")
using .Constants: Tf, ρₐ, ρₒ, ρᵢ, Cd, cᴾ, kl, Nu, α, Lat # these constants can be called inside any function
using Statistics
using KernelDensity
# change to number density * crystal volume
# turn off rise

function find_z_mean(var)
    # finds the z-mean of the array
    return vec(mean(var[1, 1, :, :], dims = 1))
end

volume = 1
aspect_ratio = 50
numSizeClasses = 200
plotstep = 1
numpdflines = 4
Rs = range(0.01, 2, numSizeClasses) .* 1e-3
Rs = exp.(range(start=log(0.01*0.001), stop=log(0.002), length=numSizeClasses))
depth = 1
Tdiff = 1e-4

ϵ = 1e-8
feltham_name = "log_space_sameConc_fast_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)

colours = cmap("CBTL1", N = Int(1.2*length(Rs))) 

cvelnum = 1
pnum = 1

feltham_filename = "1D_fields" * feltham_name * ".jld2"
end_name = feltham_name
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
if crystal_size_collision_redistribution == 2
    end_name = end_name * "_redistribute"
end


end_filename = "1D_fields" * end_name * ".jld2"
# Create a dictionary for n1 to n100 dynamically
n_dict = Dict(Symbol("n$n") => FieldTimeSeries(end_filename, "n$n") for n in 1:numSizeClasses)
# Merge with fixed fields and convert to NamedTuple
time_series = (; 
    w = FieldTimeSeries(end_filename, "w"),
    u = FieldTimeSeries(end_filename, "u"),
    v = FieldTimeSeries(end_filename, "v"),
    T = FieldTimeSeries(end_filename, "T"),
    S = FieldTimeSeries(end_filename, "S"),
    n_dict...   # Expand dictionary into NamedTuple fields
)

times = time_series.w.times
#times = times/60
#times = times/(3600*24)

# get array of data
nmatrix = zeros(numSizeClasses, length(times))
for n in range(1, stop=numSizeClasses, step=plotstep)
    nmatrix[n,:] = find_z_mean(getproperty(time_series, Symbol("n$n")))
end
nmatrix = clamp.(nmatrix, 1e-1, Inf)

########################## get the second array of data
cvelnum = 3
pnum = 2
numSizeClasses2 = 200

#feltham_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses) * "_efficiency_scattered_interpolation"

end_name1 = end_name

end_name = feltham_name
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
if crystal_size_collision_redistribution == 2
    end_name = end_name * "_redistribute"
end

end_filename = "1D_fields" * end_name * ".jld2"
# Create a dictionary for n1 to n100 dynamically
n_dict = Dict(Symbol("n$n") => FieldTimeSeries(end_filename, "n$n") for n in 1:numSizeClasses2)
# Merge with fixed fields and convert to NamedTuple
time_series2 = (; 
    w = FieldTimeSeries(end_filename, "w"),
    u = FieldTimeSeries(end_filename, "u"),
    v = FieldTimeSeries(end_filename, "v"),
    T = FieldTimeSeries(end_filename, "T"),
    S = FieldTimeSeries(end_filename, "S"),
    n_dict...   # Expand dictionary into NamedTuple fields
)

times2 = time_series2.w.times
#times2 = times2/60
#times = times/(3600*24)

nmatrix2 = zeros(numSizeClasses2, length(times2))
for n in range(1, stop=numSizeClasses2, step=plotstep)
    nmatrix2[n,:] = find_z_mean(getproperty(time_series2, Symbol("n$n")))
end
nmatrix2 = clamp.(nmatrix2, 1e-1, Inf)


# Set tick positions and labels for every 20 bars
tick_pos = [1, 40, 80, 120, 160, 200]
tick_pos = collect(range(0, stop = numSizeClasses2, step = 20))
tick_pos[1] = 1
tick_labels = string.(round.(Rs[tick_pos]*1000, digits = 2))

y_tick_pos = [1e10, 1e8, 1e6, 1e4, 1e2, 1e0] #, 1e-1]
y_tick_labels = ["10¹⁰", "10⁸", "10⁶", "10⁴", "10²", "10⁰"] #, "10⁻¹"]

labelsize = 24
ticklabelsize = 20

fig = Figure(size = (1350, 850))
ax = Axis(fig[1, 1]; xticks=(tick_pos, tick_labels), yticks=(y_tick_pos, y_tick_labels), yscale=log10, ylabel = "nᵢ", xlabel = "r(mm)", limits = ((0, numSizeClasses + 1), (1e-1, 2*maximum(nmatrix))))
bars = barplot!(ax, 1:numSizeClasses, nmatrix[:, 1], color = :gray85, fillto =0.1, gap = 0, strokecolor = :black, strokewidth = 1)

display(fig)

record(fig,  "1D_fields" * end_name * ".mp4", 1:size(nmatrix, 2)) do frame
    empty!(ax)
    ax.yscale = log10
    ax.xticks = (tick_pos, tick_labels)
    ax.yticks = (y_tick_pos, y_tick_labels)
    ax.limits = ((0, numSizeClasses + 1), (1, 2*maximum(nmatrix)))  
    #ax.limits = ((0, numSizeClasses + 1), (1e-1, 2*maximum(nmatrix)))  
    ax.ylabel = "nᵢ"
    ax.xlabel = "r (mm)" 
    ax.title = "t = " * string(round(times[frame])) * " seconds"
    ax.xlabelsize = 24
    ax.titlesize = 30
    ax.ylabelsize = 24
    ax.xticklabelsize = 20
    ax.yticklabelsize = 20
    barplot!(ax, 1:numSizeClasses, nmatrix[:, frame], color = :gray, fillto=0.1, gap = 0, strokecolor = :black, strokewidth = 1)
end

# now plot two bars on top of each other

fig = Figure(size = (1350, 850))
ax = Axis(fig[1, 1]; xticks=(tick_pos, tick_labels), yticks=(y_tick_pos, y_tick_labels), yscale=log10, ylabel = "nᵢ", xlabel = "r(mm)", limits = ((0, numSizeClasses + 1), (1e-1, 2*maximum(nmatrix))))
bars = barplot!(ax, 1:numSizeClasses, nmatrix[:, 1], color = :gray, fillto =0.1, gap = 0, strokecolor = :black, strokewidth = 1)

display(fig)

record(fig,  "1D_fields" * end_name * end_name1 * "both.mp4", 1:size(nmatrix, 2)) do frame
    empty!(ax)
    ax.yscale = log10
    ax.xticks = (tick_pos, tick_labels)
    ax.yticks = (y_tick_pos, y_tick_labels)
    ax.limits = ((0, numSizeClasses + 1), (1, 2*maximum(nmatrix)))  
    #ax.limits = ((0, numSizeClasses + 1), (1e-1, 2*maximum(nmatrix)))  
    ax.ylabel = "nᵢ"
    ax.xlabel = "r (mm)" 
    ax.title = "t = " * string(round(times[frame])) * " seconds"
    ax.xlabelsize = 24
    ax.titlesize = 30
    ax.ylabelsize = 24
    ax.xticklabelsize = 20
    ax.yticklabelsize = 20
    barplot!(ax, 1:numSizeClasses, nmatrix[:, frame], color = :gray, fillto=0.1, gap = 0, strokecolor = :black, strokewidth = 1)

    # match the times - they are matched!
    #time = times[frame]

    barplot!(ax, 1:numSizeClasses, nmatrix2[:, frame], fillto=0.1, gap = 0, color = color=RGBAf(1.0, 0.75, 0.8, 0.6), strokecolor = :red, strokewidth = 1, transparency=true)
end

c_total = zeros(length(times))
for n in range(1, stop=numSizeClasses, step=1)
    c_total = c_total .+ find_C_const(n, Rs, aspect_ratio, 1) * find_z_mean(getproperty(time_series, Symbol("n$n")))
end