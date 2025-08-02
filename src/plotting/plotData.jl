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
#Rs = exp.(range(start=log(0.01*0.001), stop=log(0.002), length=numSizeClasses))
aspect_ratio = 50
Volume = 1

#ns_to_plot = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#ns_to_plot = [1, 2, 5, 10, 50, 100, 150, 200]
#ns_to_plot = [1, 27, 62, 88, 148, 174, 189, 200]

ns_to_plot = [1, 2, 5, 10, 50, 100, 200]
#ns_to_plot = [1, 27, 62, 88, 148, 174, 200]
ns_for_color = [1, 27, 62, 88, 148, 174, 200]

#ns_to_plot = [1, 5, 40, 80, 120, 160, 200]
numpdflines = 4
colours = cmap("CBTL1", N = Int(1.2*length(Rs))) 
colours2 = cmap("CBL1", N = Int(numpdflines)+1) 


#ts1 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj_r_av.jld2", 200)
#ts2 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj_r_av_v2.jld2", 200)
#ts3 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj_r_av_v3.jld2", 200)


#ts1_times = ts1.w.times
#ts2_times = ts1_times[end] .+ ts2.w.times
#ts3_times = ts2_times[end] .+ ts3.w.times

#1D_fieldsfast_same_n_just_collisions_epsilon0.0001_200_redistribute.jld2

#ts_1 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200.jld2", numSizeClasses)
#ts_1 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj.jld2", numSizeClasses)
#ts_1 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon0.01_200.jld2", numSizeClasses)

ts_1 = load_plot_data("1D_fieldsnew_v_same_n_just_collisions_epsilon1.0e-8_200.jld2", numSizeClasses)

#ts_1 = load_plot_data("1D_fieldsnew_v_same_n_just_collisions_epsilon0.01_200.jld2", numSizeClasses)

#ts_1 = load_plot_data("1D_fieldslog_space_sameConc_fast_same_C_just_collisions_epsilon0.01_200.jld2", numSizeClasses)
#ts_1 = load_plot_data("1D_fieldslog_space_sameConc_fast_same_C_just_collisions_epsilon1.0e-8_200.jld2", numSizeClasses)

#ts_1 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_new_spherical_sum_nj.jld2", numSizeClasses) 
ts1_times = ts_1.w.times

#ts_1 = ts_nonsum
#ts1_times = ts_nonsum_times

#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon0.01_200_sum_nj.jld2", numSizeClasses) #load_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_r_av.jld2", numSizeClasses)

#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_new_spherical_redistribute.jld2", numSizeClasses)
#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj.jld2", numSizeClasses)
#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_sum_nj_redistribute.jld2", numSizeClasses)
#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon0.01_200_sum_nj.jld2", numSizeClasses)
#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon0.01_200_new_spherical.jld2", numSizeClasses) 
#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon0.01_200_new_spherical_sum_nj.jld2", numSizeClasses) 
#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon0.01_200_new_spherical_sum_nj_redistribute.jld2", numSizeClasses) 
#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_new_spherical_sum_nj_redistribute.jld2", numSizeClasses) 

ts_nonsum = load_plot_data("1D_fieldsnew_v_same_n_just_collisions_epsilon1.0e-8_200_new_spherical_sum_nj.jld2", numSizeClasses)

#ts_nonsum = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon0.01_200.jld2", numSizeClasses)

ts_nonsum_times = ts_nonsum.w.times

#ts_2 = load_plot_data("1D_fieldsfast_same_n_just_collisions_epsilon0.01_200_new_spherical_sum_nj.jld2", numSizeClasses) #load_data("1D_fieldsfast_same_n_just_collisions_epsilon1.0e-8_200_r_av.jld2", numSizeClasses)
#ts_2_times = ts_2.w.times

tmax =  2500#maximum(ts_nonsum_times)
#tmax = 700
#tmax =  maximum(ts1_times)
tmin = 0

function plot_n_C_pdf(ls, time_series_data, times, add_nC_label, add_pdf_label, end_time_series_data)
    for (indx, n) in enumerate(ns_to_plot) #range(1, stop=numSizeClasses, step=plotstep)
        if add_nC_label
            if ls == 1
                rval = round(Rs[n]*1000, digits = 2)
                lines!(ax_n, times, find_z_mean(getproperty(time_series_data, Symbol("n$n"))), linestyle=:solid, color = colours[ns_for_color[indx]], label=L"%$rval")
                #lines!(ax_n, times, find_z_mean(getproperty(time_series_data, Symbol("n$n"))), linestyle=:solid, color = colours[n], label=L"n_{%$n}")
                lines!(ax_C, times, find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(time_series_data, Symbol("n$n"))),
                        linestyle = :solid, color = colours[ns_for_color[indx]])
            end
        else
            if ls == 1
                lines!(ax_C, times, find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(time_series_data, Symbol("n$n"))),
                        linestyle = :solid, color = colours[ns_for_color[indx]])
                lines!(ax_n, times, find_z_mean(getproperty(time_series_data, Symbol("n$n"))), linestyle = :solid, color = colours[ns_for_color[indx]])
            elseif ls == 2
                lines!(ax_C, times, find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(time_series_data, Symbol("n$n"))),
                        linestyle = :dash, color = colours[ns_for_color[indx]])
                lines!(ax_n, times, find_z_mean(getproperty(time_series_data, Symbol("n$n"))), linestyle = :dash, color = colours[ns_for_color[indx]])
            else
                lines!(ax_C, times, find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(time_series_data, Symbol("n$n"))),
                        linestyle = :dot, color = colours[ns_for_color[indx]])
                lines!(ax_n, times, find_z_mean(getproperty(time_series_data, Symbol("n$n"))), linestyle = :dot, color = colours[ns_for_color[indx]])
            end
        end
    end

    # find the times when n2 has reduced by 1/10, 1/100, 1/1000, 1/10^4
    #n1max = maximum(find_z_mean(end_time_series_data.n1))
    n1max = 1.298846543545204e13
    fractions = [0.1, 0.9, 0.99, 0.999]
    #for (i, n1val) in enumerate(exp10.(range(start=log10(n1max/10), stop=log10(n1max/1.001), length=numpdflines)))
    for (i, n1val) in enumerate(fractions*n1max)
        #exp = floor(Int, log10(abs(n1val)))
        #coeff = round(n1val / 10^exp, digits = 1)
        #label_str = L"n_1 = %$coeff \times 10^{%$exp}"
        frac = fractions[i]
        #label_str = L"n_1 = %$frac n_{max}"
        label_str = L"α = %$frac"
        
        tindx = argmin(abs.(find_z_mean(time_series_data.n1) .- n1val))
        println(tindx)
        if ls == 1
            scatter!(ax_srn, times[tindx], 0.5, marker = :cross, color = colours2[i])
        else
            scatter!(ax_srn, times[tindx], 0.25,  marker = :circle, color = colours2[i])
        end
        if tindx == length(times) 
        elseif tindx == 1 #change to == 1
        else

            barvals = generatePDF(tindx, time_series_data, Rs, numSizeClasses)
            if ls == 1
                if add_pdf_label
                    #lines!(ax_rn, pdfsol.x*1000, pdfsol.density, linestyle=:solid,  color = colours[round(Int, numSizeClasses/numpdflines*i)], label = label_str)
                    lines!(ax_rn, Rs[2:end]*1000, barvals/(sum(barvals)*(Rs[1])), linestyle=:solid,  color = colours2[i], label = label_str)
                else
                    #lines!(ax_rn, pdfsol.x*1000, pdfsol.density, linestyle=:solid,  color = colours[round(Int, numSizeClasses/numpdflines*i)])
                    lines!(ax_rn, Rs[2:end]*1000, barvals/(sum(barvals)*(Rs[1])), linestyle=:solid,  color = colours2[i])
                end
            elseif ls == 2
                if add_pdf_label
                    #lines!(ax_rn, pdfsol.x*1000, pdfsol.density, linestyle=:dash,  color = colours[round(Int, numSizeClasses/numpdflines*i)], label = label_str)
                    lines!(ax_rn, Rs[2:end]*1000, barvals/(sum(barvals)*(Rs[1])), linestyle=:dash,  color = colours2[i], label = label_str)
                else
                    #lines!(ax_rn, pdfsol.x*1000, pdfsol.density, linestyle=:dash,  color = colours[round(Int, numSizeClasses/numpdflines*i)])
                    lines!(ax_rn, Rs[2:end]*1000, barvals/(sum(barvals)*(Rs[1])), linestyle=:dash,  color = colours2[i])
                end
            end
        end

    end
end

fig = Figure(size = (850, 350))

ax_srn = Axis(fig[5, 3];
            #ylabel ="n₁",
            xlabel = L"time when $n_1 = α n_{max}$ (seconds)",#)
            xgridvisible = false,
            ygridvisible = false,
            yticksvisible = false,
            yticklabelsvisible = false,
            #xscale = log10,
            #yscale = log10) #)
            #limits = ((Rs[1]*1000, Rs[end]*1000), (0, 2.2e3)))
            limits = ((tmin, tmax), (0, 0.75)))

ax_rn = Axis(fig[1:4, 3];
            ylabel =L"n(r)",
            xlabel = L"r $$(mm)",#)
            xgridvisible = false,
            ygridvisible = false,
            #xscale = log10,
            #yscale = log10, #)
            limits = ((Rs[1]*1000, Rs[end]*1000), (0, 8e3)))
            #limits = ((Rs[1]*1000, Rs[end]*1000), (1e-0, 1e3)))

ax_C = Axis(fig[1:4, 2];
            ylabel = L"C_i",
            xlabel = L"time $$(seconds)",
            yscale = log10,
            xgridvisible = false,
            ygridvisible = false,
            limits = ((tmin, tmax), (1e-14, 1e-2)))

ax_n = Axis(fig[1:4, 1];
            ylabel = L"n_i",
            xlabel = L"time $$(seconds)",
            yscale = log10,
            #limits = ((tmin, tmax), (1e4, 1e6)))
            xgridvisible = false,
            ygridvisible = false,
            limits = ((tmin, tmax), (1, 2e13)))

#plot_n_C_pdf(1, ts1, ts1_times, true, true, ts3)
#plot_n_C_pdf(1, ts2, ts2_times, false, true, ts3)
#plot_n_C_pdf(1, ts3, ts3_times, false, true, ts3)

plot_n_C_pdf(1, ts_nonsum, ts_nonsum_times, true, true, ts_nonsum)
plot_n_C_pdf(2, ts_1, ts1_times, false, false, ts_1)

#plot_n_C_pdf(2, ts_nonsum, ts_nonsum_times, false, false, ts_nonsum)
#plot_n_C_pdf(1, ts_2, ts_2_times, true, true, ts_1)

#axislegend(ax_rn, position=(:left, :top), framevisible = false)
#axislegend(ax_rn, position=(:right, :top), framevisible = false)
#axislegend(ax_rn, position=(:right, :bottom), framevisible = false)
#axislegend(ax_rn, position=(:left, :top), framevisible = false)
#axislegend(ax_C)
#axislegend(ax_n, "rᵢ (mm)", framevisible = false)
#axislegend(ax_n, "rᵢ (mm)", framevisible = false, nbanks = 2)

elem_1 = [LineElement(color = :red, linestyle = nothing),
          MarkerElement(color = :blue, marker = 'x', markersize = 15,
          strokecolor = :black)]

elem_2 = [LineElement(color = :black, linestyle = :dash)]

#Legend(f[1, 2],
#    [elem_1, elem_2, elem_3, elem_4, elem_5],
#    ["Line & Marker", "Poly & Line", "Line", "Marker", "Poly"],
#    patchsize = (35, 35), rowgap = 10)

scatterlines!(ax_n, [-1, -2], [1 , 1], linestyle=:solid, marker = :circle, color = :black, label= "old model")
scatterlines!(ax_n, [-1, -2], [1 , 1], linestyle=:dash, color = :black, marker = :cross, label= "new model")
Legend(fig[5, 1:2], ax_n, "rᵢ (mm)", orientation = :horizontal, nbanks = 2, titleposition = :left, framevisible = false, padding = (0.0f0, 0.0f0, 0.0f0, 0.0f0))

fac = 2/3
letters = ["(a)", "(b)", "(c)", "(d)"]
    for i = 1:4
        if i < 3
            label_a = fig[1:4, i, TopLeft()] = Label(fig, letters[i], fontsize = 24*fac, halign = :right)
            label_a.padding = (0, 40, -10, 0)
        elseif i == 3
            label_a = fig[1:4, i, TopLeft()] = Label(fig, letters[i], fontsize = 24*fac, halign = :right)
            label_a.padding = (0, 40, -10, 0)
        else
            label_b = fig[5, 3, TopLeft()] = Label(fig, letters[i], fontsize = 24*fac, halign = :right)
            label_b.padding = (0, 40, -10, 0)
        end
    end

#for (ax, label) in zip([ax_n, ax_C, ax_rn, ax_srn], ["(a)", "(b)", "(c)", "(d)"])
#    text!(
#        ax, 0, 1,
#        text = label,
#        #font = :bold,
#        align = (:left, :top),
#        offset = (6, -2),
#        space = :relative,
#        fontsize = 16
#    )
#end

fig

#save("dot_interact_eps0.01_c13.pdf", fig)
#save("dot_interact_eps0.01_c3.pdf", fig)
#save("n_sum_redistribute_eps0.01_c3.pdf", fig)
#save("dot_interact_eps0.01_c1.pdf", fig)

#save("nsumresdistribute_eps1e-8.pdf", fig)
#save("ndot_eps1e-8.pdf", fig)
save("nnewdot_eps1e-8.pdf", fig)
#save("test.pdf", fig)


#Vs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] * 1e-9
#Rs = ( (Vs * aspect_ratio) / (2 * π) ) .^(1/3) 

c_total = zeros(length(ts_nonsum_times))
for n in range(1, stop=numSizeClasses, step=1)
    c_total = c_total .+ find_C_const(n, Rs, aspect_ratio, Volume) * find_z_mean(getproperty(ts_nonsum, Symbol("n$n")))
end
c_total