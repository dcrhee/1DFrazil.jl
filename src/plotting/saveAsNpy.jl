using JLD2
using NPZ
include("loadData.jl")
using .loadData

include("plottingFunctions.jl")
using .plottingFunctions

for ϵ = [1e-3, 1e-4, 1e-5, 1e-8, 1e-2]

numSizeClasses = 200
Rs = range(0.01, 2, numSizeClasses) .* 1e-3
#Rs = exp.(range(start=log(0.01*0.001), stop=log(0.002), length=numSizeClasses))
aspect_ratio = 50
Volume = 1

for crystal_size_collision_redistribution = [2]
for cvelnum = [3]#[1, 3]
for pnum = [2]#[1, 2ß]

#end_name = "new_v_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
#end_name = "new_v_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses) * "_new_spherical_sum_nj_r_max_redistribute"
end_name = "new_v_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
#end_name = "new_v_log_space_same_C_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
new_label = ""

collision_velocity_parameterisation_num = cvelnum # 1 is old cylinder, 2 is new cylinder, 3 is new spherical
concentration_parameterisation_num = pnum # 1 is mean n, 2 is sum over nj
effective_radius = false # add in their effective radius
efficiency_radius = false
max_radius = true

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
if max_radius
    end_name = end_name * "_r_max"
end
if efficiency_radius
    end_name = end_name * "_r_av"
end
if crystal_size_collision_redistribution == 2
    end_name = end_name * "_redistribute"
end

end_filename = "1D_fields" * end_name * ".jld2"
try
    println("end_name", end_filename)
    ts = load_plot_data(end_filename, numSizeClasses)
    ts_times = ts.w.times

    ns = zeros(numSizeClasses + 1, length(ts_times))
    ns[1, :] = ts_times
    for nindx in 1:numSizeClasses
        ns[nindx+1, :] = find_z_mean(getproperty(ts, Symbol("n$nindx")))
    end

    npzwrite(end_name * ".npy", ns)

catch
end
end

end
end
end