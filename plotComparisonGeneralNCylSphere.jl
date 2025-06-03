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

function find_C_const(indx)
    # finds the z-mean of the array
    return 2*π*Rs[indx]^3/(aspect_ratio*volume)
end

#Rs = [0.01, 0.05, 0.15, 0.3, 0.4, 0.5, 0.6, 0.8, 1, 2] .* 1e-3
volume = 1
aspect_ratio = 50
numSizeClasses = 200
plotstep = 39
numpdflines = 4
Rs = range(0.01, 2, numSizeClasses) .* 1e-3
depth = 1

ϵ = 1e-2 #7.4 * 1e-6 #10^-3 # m²s⁻³
#feltham_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)#"_same_n_just_collisions_epsilon" * string(ϵ)
feltham_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses) * "_depth_" * string(depth)
#feltham_name = "fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)  * "_sum_nj"


#colours = cmap("Gouldian", N = length(Rs))
colours = cmap("CBTL1", N = Int(1.2*length(Rs))) 

# Core variables
feltham_filename = "1D_fields" * feltham_name * ".jld2"
# Create a dictionary for n1 to n100 dynamically
n_dict = Dict(Symbol("n$n") => FieldTimeSeries(feltham_filename, "n$n") for n in 1:numSizeClasses)
# Merge with fixed fields and convert to NamedTuple
feltham_time_series = (; 
    w = FieldTimeSeries(feltham_filename, "w"),
    u = FieldTimeSeries(feltham_filename, "u"),
    v = FieldTimeSeries(feltham_filename, "v"),
    T = FieldTimeSeries(feltham_filename, "T"),
    S = FieldTimeSeries(feltham_filename, "S"),
    n_dict...   # Expand dictionary into NamedTuple fields
)

feltham_times = feltham_time_series.w.times
feltham_times = feltham_times/60
#feltham_times = feltham_times/(3600*24)

cvelnum = 3
pnum = 2
#for cvelnum = [1]#[1, 2, 3]
#for pnum = [2]#[1, 2]
end_name = feltham_name #"fast_same_n_just_collisions_epsilon" * string(ϵ) * "_" * string(numSizeClasses)
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
times = times/60
#times = times/(3600*24)


tmax =  maximum(times)
tmin = 0

fig = Figure(size = (850, 850))

ax_no = Axis(fig[2, 1:2];
            ylabel = "n, mean density",
            xlabel = "time (minutes)",
            yscale = log10,
            #limits = ((tmin, 4), (0.1, 2e13)))
            limits = ((tmin, tmax), (0.1, 2e13)))

ax_C = Axis(fig[1, 3:4];
            ylabel = "C",
            xlabel = "time (minutes)",
            yscale = log10,
            limits = ((tmin, tmax), (1e-15, 2e-3)))

ax_Co = Axis(fig[1, 1:2],
            ylabel = "C, mean density",
            xlabel = "time (minutes)",
            yscale = log10,
            #limits = ((tmin, 4), (1e-15, 2e-3)))
            limits = ((tmin, tmax), (1e-15, 2e-3)))

ax_n = Axis(fig[2, 3:4];
            ylabel = "n",
            xlabel = "time (minutes)",
            yscale = log10,
            #limits = ((tmin, tmax), (1e4, 1e6)))
            limits = ((tmin, tmax), (0.1, 2e13)))

for n in range(1, stop=numSizeClasses, step=plotstep)
    lines!(ax_Co, feltham_times, find_C_const(n) * find_z_mean(getproperty(feltham_time_series, Symbol("n$n"))),
            linestyle=:solid, color=colours[n])
    lines!(ax_no, feltham_times, find_z_mean(getproperty(feltham_time_series, Symbol("n$n"))), linestyle=:solid, color = colours[n])

    lines!(ax_C, times, find_C_const(n) * find_z_mean(getproperty(time_series, Symbol("n$n"))),
            linestyle=:solid, color=colours[n], label=L"C_{%$n}")
    lines!(ax_n, times, find_z_mean(getproperty(time_series, Symbol("n$n"))), linestyle=:solid, color = colours[n], label=L"n_{%$n}")
end

axislegend(ax_C)
axislegend(ax_n)

fig

save(end_name * "_compare_nomax.png", fig)
save(end_name * "_compare_nomax.eps", fig)

# plot the figures on top of each other

# also plot the mean number and mean concentration
# plot mean radius against time (one line including smallest size class, other line excluding smallest size class)

# Efficient computation of weighted sum
mean_radius_f = sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) .* Rs[i] for i in 1:numSizeClasses)./sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) for i in 1:numSizeClasses)
mean_radius_exclude_small_f = sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) .* Rs[i] for i in 2:numSizeClasses)./sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) for i in 2:numSizeClasses)
num_g1_f = sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) for i in 2:numSizeClasses)

mean_radius = sum(find_z_mean(getfield(time_series, Symbol("n$i"))) .* Rs[i] for i in 1:numSizeClasses)./sum(find_z_mean(getfield(time_series, Symbol("n$i"))) for i in 1:numSizeClasses)
mean_radius_exclude_small = sum(find_z_mean(getfield(time_series, Symbol("n$i"))) .* Rs[i] for i in 2:numSizeClasses)./sum(find_z_mean(getfield(time_series, Symbol("n$i"))) for i in 2:numSizeClasses)
num_g1 = sum(find_z_mean(getfield(time_series, Symbol("n$i"))) for i in 2:numSizeClasses)

fig2 = Figure(size = (850, 450))

ax_ro = Axis(fig2[1, 1];
            ylabel ="r̅",
            xlabel = "time (days)",
            yscale = log10,
            limits = ((tmin, tmax), (9e-6, 1e-2)))

ax_r = Axis(fig2[1, 2];
            ylabel = "r̅",
            xlabel = "time (days)",
            yscale = log10,
            limits = ((tmin, tmax), (9e-6, 1e-2)))


lines!(ax_ro, feltham_times, mean_radius_f, linestyle=:solid,  color = colours[1], label = "all i")
lines!(ax_ro, feltham_times[num_g1_f .> 0], mean_radius_exclude_small_f[num_g1_f .> 0], linestyle=:solid,  color = colours[numSizeClasses-1], label = "i > 1")

lines!(ax_r, times, mean_radius, linestyle=:solid,  color = colours[1], label = "all i")
lines!(ax_r, times[num_g1 .> 0], mean_radius_exclude_small[num_g1 .> 0], linestyle=:solid,  color = colours[numSizeClasses-1], label = "i > 1")

axislegend(ax_ro)

fig2
save(end_name * "_compare_rs_zoom.png", fig2)
save(end_name * "_compare_rs_zoom.eps", fig2)

conc = sum(find_z_mean(getfield(feltham_time_series, Symbol("n$i"))) .* 2*π*Rs[i]^3/50 for i in 1:numSizeClasses)
conc = sum(find_z_mean(getfield(time_series, Symbol("n$i"))) .* 2*π*Rs[i]^3/50 for i in 1:numSizeClasses)

########################################### histograms ###################################################

# generate pdf
function generatePDF(tindx, chosentimeseries)
    # input: tindx = time
    expanded_data = Vector{Float64}()  # Pre-allocate an empty vector
    barvals = zeros(length(Rs)-1)
    for i in 2:numSizeClasses
        nivals = find_z_mean(getfield(chosentimeseries, Symbol("n$i")))  # Extract tracer values
        append!(expanded_data, Rs[i] * ones(round(Int, nivals[tindx])))
        barvals[i-1] = nivals[tindx]
    end
        
    # Perform Kernel Density Estimation of everything apart from size class 1
    pdf = kde(expanded_data)

    return pdf, barvals
end
# Create a dataset where each r is repeated according to its n(r)

# plot the histogram and then on top of it plot the density estimate from the pdf
fig3 = Figure(size = (850, 450))

ax_h = Axis(fig3[1, 1];
            ylabel ="r̅",
            xlabel = "time (days)",
            yscale = log10,)

ax_p = Axis(fig3[1, 2];
            ylabel = "PDF",
            xlabel = "r",
            yscale = log10,)

tindx = 50
pdf, barvals = generatePDF(tindx, time_series)
lines!(ax_p, pdf.x, pdf.density, linestyle=:solid,  color = colours[1], label = "all i")
barplot!(ax_h, Rs[2:end], barvals)

fig3

#################################################

fig = Figure(size = (1050, 450))

ax_r = Axis(fig[1, 5:6];
            ylabel ="n(r)",
            xlabel = "r (mm)")
            #yscale = log10)

ax_C = Axis(fig[1, 3:4];
            ylabel = "Cᵢ",
            xlabel = "time (minutes)",
            yscale = log10,
            limits = ((tmin, tmax), (1e-15, 2e-3)))

ax_n = Axis(fig[1, 1:2];
            ylabel = "nᵢ",
            xlabel = "time (minutes)",
            yscale = log10,
            #limits = ((tmin, tmax), (1e4, 1e6)))
            limits = ((tmin, tmax), (0.1, 2e13)))

# find the times when n2 has reduced by 1/10, 1/100, 1/1000, 1/10^4
n1max = maximum(find_z_mean(time_series.n1))
n1initial = time_series.n1[1, 1, 1, 1]

plottimes = [5.4, 5.9, 6.5, 10.0]

#for (i, n1val) in enumerate(exp10.(range(start=log10(n1initial), stop=log10(n1max/1.1), length=plotstep)))
for (i, n1val) in enumerate(exp10.(range(start=log10(n1max/10), stop=log10(n1max/1.001), length=numpdflines)))
    tindx = argmin(abs.(find_z_mean(time_series.n1) .- n1val))
    tindx = argmin(abs.(plottimes[i] .- feltham_times))
    try
        pdf, barvals = generatePDF(tindx, feltham_time_series)
        lines!(ax_r, pdf.x*1000, pdf.density, linestyle=:dash,  color = colours[round(Int, numSizeClasses/numpdflines*i)])
    catch
    end
    pdf, barvals = generatePDF(tindx, time_series)
    lines!(ax_r, pdf.x*1000, pdf.density, linestyle=:solid,  color = colours[round(Int, numSizeClasses/numpdflines*i)], label = "t = " * string(feltham_times[tindx]) * " minutes")
end

#lines!(ax_r, feltham_times, mean_radius_f, linestyle=:dash,  color = colours[1])
#lines!(ax_r, feltham_times[num_g1_f .> 0], mean_radius_exclude_small_f[num_g1_f .> 0], linestyle=:dash,  color = colours[9])

#lines!(ax_r, times, mean_radius, linestyle=:solid,  color = colours[1], label = "all i")
#lines!(ax_r, times[num_g1 .> 0], mean_radius_exclude_small[num_g1 .> 0], linestyle=:solid,  color = colours[9], label = "i > 1")

axislegend(ax_r, position=(:right, :bottom))

for n in range(1, stop=numSizeClasses, step=plotstep)
    lines!(ax_C, feltham_times, find_C_const(n) * find_z_mean(getproperty(feltham_time_series, Symbol("n$n"))),
            linestyle=:dash, color=colours[n])
   lines!(ax_n, feltham_times, find_z_mean(getproperty(feltham_time_series, Symbol("n$n"))), linestyle=:dash, color = colours[n])

    lines!(ax_C, times, find_C_const(n) * find_z_mean(getproperty(time_series, Symbol("n$n"))),
            linestyle=:solid, color=colours[n], label=L"C_{%$n}")
    lines!(ax_n, times, find_z_mean(getproperty(time_series, Symbol("n$n"))), linestyle=:solid, color = colours[n], label=L"n_{%$n}")
end

axislegend(ax_C)
axislegend(ax_n)

fig

#save(end_name * "_4_compare_zoom_max.png", fig)
save(end_name * "_compare_nomax_ontop.png", fig)
save(end_name * "_compare_nomax_ontop.eps", fig)

fig = Figure(size = (1050, 450))

ax_r = Axis(fig[1, 5:6];
            ylabel ="n(r)",
            xlabel = "r (mm)")
            #yscale = log10)

ax_C = Axis(fig[1, 3:4];
            ylabel = "Cᵢ",
            xlabel = "time (minutes)",
            yscale = log10,
            limits = ((tmin, tmax), (1e-15, 2e-3)))

ax_n = Axis(fig[1, 1:2];
            ylabel = "nᵢ",
            xlabel = "time (minutes)",
            yscale = log10,
            #limits = ((tmin, tmax), (1e4, 1e6)))
            limits = ((tmin, tmax), (0.1, 2e13)))

# find the times when n2 has reduced by 1/10, 1/100, 1/1000, 1/10^4
n1max = maximum(find_z_mean(time_series.n1))
n1initial = time_series.n1[1, 1, 1, 1]


for (i, n1val) in enumerate(exp10.(range(start=log10(n1max/10), stop=log10(n1max/1.001), length=numpdflines)))
    exp = floor(Int, log10(abs(n1val)))
    coeff = round(n1val / 10^exp, digits = 1)
    label_str = L"n_1 = %$coeff \times 10^{%$exp}"

    tindx = argmin(abs.(find_z_mean(feltham_time_series.n1) .- n1val))
    pdf, barvals = generatePDF(tindx, feltham_time_series)
    lines!(ax_r, pdf.x*1000, pdf.density, linestyle=:dash,  color = colours[round(Int, numSizeClasses/numpdflines*i)])

    tindx = argmin(abs.(find_z_mean(time_series.n1) .- n1val))
    pdf, barvals = generatePDF(tindx, time_series)
    lines!(ax_r, pdf.x*1000, pdf.density, linestyle=:solid,  color = colours[round(Int, numSizeClasses/numpdflines*i)], label = label_str)
end
for n in range(1, stop=numSizeClasses, step=plotstep)
    lines!(ax_C, feltham_times, find_C_const(n) * find_z_mean(getproperty(feltham_time_series, Symbol("n$n"))),
            linestyle=:dash, color=colours[n])
   lines!(ax_n, feltham_times, find_z_mean(getproperty(feltham_time_series, Symbol("n$n"))), linestyle=:dash, color = colours[n])

    lines!(ax_C, times, find_C_const(n) * find_z_mean(getproperty(time_series, Symbol("n$n"))),
            linestyle=:solid, color=colours[n], label=L"C_{%$n}")
    lines!(ax_n, times, find_z_mean(getproperty(time_series, Symbol("n$n"))), linestyle=:solid, color = colours[n], label=L"n_{%$n}")
end
axislegend(ax_r, position=(:right, :bottom))
axislegend(ax_C)
axislegend(ax_n)

fig

#save(end_name * "_4_compare_zoom_max.png", fig)
save(end_name * "_compare_nomax_ontop_nr.png", fig)
save(end_name * "_compare_nomax_ontop_nr.eps", fig)


fig = Figure(size = (850, 850))

ax_rn = Axis(fig[2, 1];
            ylabel ="n(r)",
            xlabel = "r (mm)")
            #yscale = log10)

ax_rt = Axis(fig[2, 2];
            ylabel ="n(r)",
            xlabel = "r (mm)")
            #yscale = log10)

ax_C = Axis(fig[1, 2];
            ylabel = "Cᵢ",
            xlabel = "time (minutes)",
            yscale = log10,
            limits = ((tmin, tmax), (1e-15, 2e-3)))

ax_n = Axis(fig[1, 1];
            ylabel = "nᵢ",
            xlabel = "time (minutes)",
            yscale = log10,
            #limits = ((tmin, tmax), (1e4, 1e6)))
            limits = ((tmin, tmax), (0.1, 2e13)))

# find the times when n2 has reduced by 1/10, 1/100, 1/1000, 1/10^4
n1max = maximum(find_z_mean(time_series.n1))
n1initial = time_series.n1[1, 1, 1, 1]

#for (i, n1val) in enumerate(exp10.(range(start=log10(n1initial), stop=log10(n1max/1.1), length=plotstep)))
for (i, n1val) in enumerate(exp10.(range(start=log10(n1max/10), stop=log10(n1max/1.001), length=numpdflines)))
    tindx = argmin(abs.(find_z_mean(time_series.n1) .- n1val))
    try
        pdf, barvals = generatePDF(tindx, feltham_time_series)
        lines!(ax_rt, pdf.x*1000, pdf.density, linestyle=:dash,  color = colours[round(Int, numSizeClasses/numpdflines*i)])
    catch
    end

    pdf, barvals = generatePDF(tindx, time_series)
    lines!(ax_rt, pdf.x*1000, pdf.density, linestyle=:solid,  color = colours[round(Int, numSizeClasses/numpdflines*i)], label = "t = " * string(feltham_times[tindx]) * " minutes")
end

for (i, n1val) in enumerate(exp10.(range(start=log10(n1max/10), stop=log10(n1max/1.001), length=numpdflines)))
    exp = floor(Int, log10(abs(n1val)))
    coeff = round(n1val / 10^exp, digits = 1)
    label_str = L"n_1 = %$coeff \times 10^{%$exp}"

    tindx = argmin(abs.(find_z_mean(feltham_time_series.n1) .- n1val))
    pdf, barvals = generatePDF(tindx, feltham_time_series)
    lines!(ax_rn, pdf.x*1000, pdf.density, linestyle=:dash,  color = colours[round(Int, numSizeClasses/numpdflines*i)])

    tindx = argmin(abs.(find_z_mean(time_series.n1) .- n1val))
    pdf, barvals = generatePDF(tindx, time_series)
    lines!(ax_rn, pdf.x*1000, pdf.density, linestyle=:solid,  color = colours[round(Int, numSizeClasses/numpdflines*i)], label = label_str)
end

#lines!(ax_r, feltham_times, mean_radius_f, linestyle=:dash,  color = colours[1])
#lines!(ax_r, feltham_times[num_g1_f .> 0], mean_radius_exclude_small_f[num_g1_f .> 0], linestyle=:dash,  color = colours[9])

#lines!(ax_r, times, mean_radius, linestyle=:solid,  color = colours[1], label = "all i")
#lines!(ax_r, times[num_g1 .> 0], mean_radius_exclude_small[num_g1 .> 0], linestyle=:solid,  color = colours[9], label = "i > 1")

axislegend(ax_rt, position=(:right, :bottom))
axislegend(ax_rn, position=(:right, :bottom))

for n in range(1, stop=numSizeClasses, step=plotstep)
    lines!(ax_C, feltham_times, find_C_const(n) * find_z_mean(getproperty(feltham_time_series, Symbol("n$n"))),
            linestyle=:dash, color=colours[n])
    lines!(ax_n, feltham_times, find_z_mean(getproperty(feltham_time_series, Symbol("n$n"))), linestyle=:dash, color = colours[n])

    lines!(ax_C, times, find_C_const(n) * find_z_mean(getproperty(time_series, Symbol("n$n"))),
            linestyle=:solid, color=colours[n], label=L"C_{%$n}")
    lines!(ax_n, times, find_z_mean(getproperty(time_series, Symbol("n$n"))), linestyle=:solid, color = colours[n], label=L"n_{%$n}")
end

axislegend(ax_C)
axislegend(ax_n)

fig

#save(end_name * "_4_compare_zoom_max.png", fig)
save(end_name * "_compare_nomax_4ontop.png", fig)
save(end_name * "_compare_nomax_4ontop.eps", fig)
