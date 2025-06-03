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

numSizeClasses = 200
Rs = range(0.01, 2, numSizeClasses) .* 1e-3
aspect_ratio = 50
Volume = 1

#ns_to_plot = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

ns_to_plot = [1, 40, 80, 120, 160, 200]
numpdflines = 4
colours = cmap("CBTL1", N = Int(1.2*length(Rs))) 

#ts1 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj_r_av.jld2", 200)
#ts2 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj_r_av_v2.jld2", 200)
#ts3 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj_r_av_v3.jld2", 200)


#ts1_times = ts1.w.times
#ts2_times = ts1_times[end] .+ ts2.w.times
#ts3_times = ts2_times[end] .+ ts3.w.times

#tssphere = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon0.01_200_new_spherical_sum_nj_r_av.jld2", numSizeClasses)
#sphere_times = tssphere.w.times

ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj_redistribute.jld2", numSizeClasses) #load_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_r_av.jld2", numSizeClasses)
ts_nonsum_times = ts_nonsum.w.times

#tmax =  maximum(vcat(ts1_times, ts2_times, ts3_times))
tmax =  maximum(ts_nonsum_times)
tmin = 0

function plot_n_C_pdf(ls, time_series_data, times, add_nC_label, add_pdf_label, end_time_series_data)
    for n in ns_to_plot #range(1, stop=numSizeClasses, step=plotstep)
        if add_nC_label
            if ls == 1
                lines!(ax_n, times, find_z_mean(getproperty(time_series_data, Symbol("n$n"))), linestyle=:solid, color = colours[n], label=L"n_{%$n}")
                lines!(ax_C, times, find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(time_series_data, Symbol("n$n"))),
                        linestyle = :solid, color=colours[n])
            end
        else
            if ls == 1
                lines!(ax_C, times, find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(time_series_data, Symbol("n$n"))),
                        linestyle = :solid, color=colours[n])
                lines!(ax_n, times, find_z_mean(getproperty(time_series_data, Symbol("n$n"))), linestyle = :solid, color = colours[n])
            elseif ls == 2
                lines!(ax_C, times, find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(time_series_data, Symbol("n$n"))),
                        linestyle = :dash, color=colours[n])
                lines!(ax_n, times, find_z_mean(getproperty(time_series_data, Symbol("n$n"))), linestyle = :dash, color = colours[n])
            else
                lines!(ax_C, times, find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(time_series_data, Symbol("n$n"))),
                        linestyle = :dot, color=colours[n])
                lines!(ax_n, times, find_z_mean(getproperty(time_series_data, Symbol("n$n"))), linestyle = :dot, color = colours[n])
            end
        end
    end

    # find the times when n2 has reduced by 1/10, 1/100, 1/1000, 1/10^4
    n1max = maximum(find_z_mean(end_time_series_data.n1))
    
    for (i, n1val) in enumerate(exp10.(range(start=log10(n1max/10), stop=log10(n1max/1.001), length=numpdflines)))
        exp = floor(Int, log10(abs(n1val)))
        coeff = round(n1val / 10^exp, digits = 1)
        label_str = L"n_1 = %$coeff \times 10^{%$exp}"

        tindx = argmin(abs.(find_z_mean(time_series_data.n1) .- n1val))
        if tindx == length(times) 
        elseif tindx  == 1
        else

            pdfsol, barvals = generatePDF(tindx, time_series_data, Rs, numSizeClasses)
            if ls == 1
                if add_pdf_label
                    lines!(ax_rn, pdfsol.x*1000, pdfsol.density, linestyle=:solid,  color = colours[round(Int, numSizeClasses/numpdflines*i)], label = label_str)
                else
                    lines!(ax_rn, pdfsol.x*1000, pdfsol.density, linestyle=:solid,  color = colours[round(Int, numSizeClasses/numpdflines*i)])
                end
            elseif ls == 2
                if add_pdf_label
                    lines!(ax_rn, pdfsol.x*1000, pdfsol.density, linestyle=:dash,  color = colours[round(Int, numSizeClasses/numpdflines*i)], label = label_str)
                else
                    lines!(ax_rn, pdfsol.x*1000, pdfsol.density, linestyle=:dash,  color = colours[round(Int, numSizeClasses/numpdflines*i)])
                end
            end
        end

    end
end

fig = Figure(size = (850, 350))

ax_rn = Axis(fig[1, 3];
            ylabel ="n(r)",
            xlabel = "r (mm)")
            #yscale = log10)

ax_C = Axis(fig[1, 2];
            ylabel = "Cᵢ",
            xlabel = "time (seconds)",
            yscale = log10,
            limits = ((tmin, tmax), (1e-15, 1e-2)))

ax_n = Axis(fig[1, 1];
            ylabel = "nᵢ",
            xlabel = "time (seconds)",
            yscale = log10,
            #limits = ((tmin, tmax), (1e4, 1e6)))
            limits = ((tmin, tmax), (1, 2e13)))

#plot_n_C_pdf(1, ts1, ts1_times, true, true, ts3)
#plot_n_C_pdf(1, ts2, ts2_times, false, true, ts3)
#plot_n_C_pdf(1, ts3, ts3_times, false, true, ts3)

plot_n_C_pdf(1, ts_nonsum, ts_nonsum_times, true, true, ts_nonsum)
#plot_n_C_pdf(2, tssphere, sphere_times, false, false, tssphere)



axislegend(ax_rn, position=(:left, :top))
#axislegend(ax_C)
axislegend(ax_n)


fig

#save(end_name * "_4_compare_zoom_max.png", fig)
#save("r_rav_eps1e-8.png", fig)
#save("compare_rsum_rav_eps1e-8.png", fig)
#save(end_name * "_compare_nomax_3ontop.eps", fig)

#Vs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] * 1e-9
#Rs = ( (Vs * aspect_ratio) / (2 * π) ) .^(1/3) 

c_total = zeros(length(ts_nonsum_times))
for n in range(1, stop=numSizeClasses, step=1)
    c_total = c_total .+ find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(ts_nonsum, Symbol("n$n")))
end
