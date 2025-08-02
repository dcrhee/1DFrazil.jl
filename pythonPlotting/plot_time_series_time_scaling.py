#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 26 11:25:04 2025

@author: cotton
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

def find_time_scale_factor(ns1, ns2):
    fractions_to_match = np.array([1e-3, 1e-2, 0.1, 0.9, 0.99, 0.999])
    n1_val_to_match = n1max * fractions_to_match
    
    # create interpolated curve from n1 > n1initial
    n1_initial = 1e9
    
    # create interpolated curves
    intepolated_ns1 = interpolate.CubicSpline(ns1[0, np.sum(ns1[1, :] <= n1_initial):], np.log10(ns1[1, np.sum(ns1[1, :] <= n1_initial):]))
    intepolated_ns2 = interpolate.CubicSpline(ns2[0, np.sum(ns2[1, :] <= n1_initial):], np.log10(ns2[1, np.sum(ns2[1, :] <= n1_initial):]))
    numpoints = 1000
    
    n1_times = np.zeros(len(n1_val_to_match))
    n2_times = np.zeros(len(n1_val_to_match))
    for indx, n1_val in enumerate(n1_val_to_match):
        
        
        t1_start_indx = np.sum(ns1[1, :] < n1_val)-1
        t1_end_indx= t1_start_indx + 1
        
        t2_start_indx = np.sum(ns2[1, :] < n1_val)-1
        t2_end_indx= t2_start_indx + 1
        
        n1_times_to_interp = np.linspace(ns1[0, t1_start_indx], ns1[0, t1_end_indx], numpoints)
        n2_times_to_interp = np.linspace(ns2[0, t2_start_indx], ns2[0, t2_end_indx], numpoints)
        
        n1_interpolate_n1vals = intepolated_ns1(n1_times_to_interp)

        n1_times[indx] = n1_times_to_interp[np.argmin(np.abs(n1_interpolate_n1vals - np.log10(n1_val)))]
        
        n2_interpolate_n1vals = intepolated_ns2(n2_times_to_interp)

        n2_times[indx] = n2_times_to_interp[np.argmin(np.abs(n2_interpolate_n1vals - np.log10(n1_val)))]
        
        
    #axm["n"].semilogy(np.linspace(ns1[0, np.sum(ns1[1, :] <= n1_initial)], ns1[0, -1], 10000)*time_scale_factor, 10**intepolated_ns1(np.linspace(ns1[0, np.sum(ns1[1, :] <= n1_initial)], ns1[0, -1], 10000)), '.', color = "y")    
    #axm["n"].semilogy(n1_times*time_scale_factor, n1_val_to_match, 'o', color = "y")
    #axm["n"].semilogy(n2_times, n1_val_to_match, '.', color = "b")
    return n1_times, n2_times


def plot_data(axm, SD, num, typelabel, time_scale_factor):
    
    if num == 1:          
        
        for indx, n in enumerate(SD_to_plot):
            axm["n"].semilogy(SD[0, :]*time_scale_factor, SD[n, :], color = colours[indx], label = r"$n_{" + str(n) + "}$")
            axm["c"].semilogy(SD[0, :]*time_scale_factor, cfact[n-1]*SD[n, :], color = colours[indx], label = r"$C_{" + str(n) + "}$")
     
        
        for i, n1val in enumerate(fractions*n1max):
            frac = fractions[i]
            label_str = r"$\alpha = " + str(frac) + "$"
            
            tindx = np.argmin(np.abs(SD[1, :]- n1val))
            if initial_rad == "same_n":
                axm["p"].plot(Rs[1:]*1000, SD[2:, tindx]/np.sum(SD[2:, tindx]), color = colours2[(i+1)], ls = "-", label = label_str)
            elif initial_rad == "log_space_same_C":
                axm["p"].semilogx(Rs[1:]*1000, cfact[1:]*SD[2:, tindx]/np.sum(cfact[1:]*SD[2:, tindx]), color = colours2[(i+1)], ls = "-", label = label_str)
            else:
                axm["p"].plot(Rs[1:]*1000, cfact[1:]*SD[2:, tindx]/np.sum(cfact[1:]*SD[2:, tindx]), color = colours2[(i+1)], ls = "-", label = label_str)
            axm["t"].plot(SD[0, tindx], [0.5], 'x', color = colours2[(i+1)])
        print("SDR", tindx)
        #axm["n"].semilogy([0, 0], [-1, -1], color = "k", ls = "-", label = typelabel)
        #axm["c"].semilogy([0, 0], [-1, -1], color = "k", ls = "-", label = typelabel)
        if initial_rad == "same_n":
            #axm["p"].plot(Rs[1:]*1000-10, SD[2:, tindx]/np.sum(SD[2:, tindx]), "-", color = "k", label = typelabel)
            axm["t"].plot([-100], [0.5], 'x', color = "k", label = typelabel)
            axm["t"].plot([-100], [0.5], '-', color = "k", label = typelabel)
        else:
            #axm["p"].plot(Rs[1:]*1000-10, cfact[1:]*SD[2:, tindx]/np.sum(cfact[1:]*SD[2:, tindx]), "-", color = "k", label = typelabel)
            axm["t"].plot([-100], [0.5], 'x', color = "k", label = typelabel)
            axm["t"].plot([-100], [0.5], '-', color = "k", label = typelabel)
        

    else:
        #axm["n"].semilogy([0, 0], [-1, -1], ls = "--", color = "k", label = typelabel)
        #axm["c"].semilogy([0, 0], [-1, -1], ls = "--", color = "k", label = typelabel)
        for indx, n in enumerate(SD_to_plot):
            axm["n"].semilogy(SD[0, :], SD[n, :], color = colours[indx], ls = "--")
            axm["c"].semilogy(SD[0, :], cfact[n-1]*SD[n, :], color = colours[indx], ls = "--")
    
        
        for i, n1val in enumerate(fractions*n1max):
            frac = fractions[i]
            label_str = r"$\alpha = " + str(frac) + "$"
            
            tindx = np.argmin(np.abs(SD[1, :]- n1val))
            if initial_rad == "same_n":
                axm["p"].plot(Rs[1:]*1000, SD[2:, tindx]/np.sum(SD[2:, tindx]), color = colours2[(i+1)], ls = "--")
            elif initial_rad == "log_space_same_C":
                 axm["p"].semilogx(Rs[1:]*1000, cfact[1:]*SD[2:, tindx]/np.sum(cfact[1:]*SD[2:, tindx]), color = colours2[(i+1)], ls = "--")                
            else:
                axm["p"].plot(Rs[1:]*1000, cfact[1:]*SD[2:, tindx]/np.sum(cfact[1:]*SD[2:, tindx]), color = colours2[(i+1)], ls = "--")
            axm["t"].plot(SD[0, tindx], [-0.5], 'o', color = colours2[(i+1)])
        print("SD", tindx)
        
        if initial_rad == "same_n":
            #axm["p"].plot(Rs[1:]*1000-10, SD[2:, tindx]/np.sum(SD[2:, tindx]), "--", color = "k", label = typelabel)
            axm["t"].plot([-100], [0.5], 'o', color = "k", label = typelabel)
            axm["t"].plot([-100], [0.5], '--', color = "k", label = typelabel)
        else:
            #axm["p"].plot(Rs[1:]*1000-10, cfact[1:]*SD[2:, tindx]/np.sum(cfact[1:]*SD[2:, tindx]), "--", color = "k", label = typelabel)
            axm["t"].plot([-100], [-0.5], 'o', color = "k", label = typelabel)
            axm["t"].plot([-100], [-0.5], '--', color = "k", label = typelabel)

    return

fig, axm = plt.subplot_mosaic(mosaic= """
ncp
ncp
ncp
ncp
ncp
ncp
nct
nct
""", layout="tight", figsize=(10, 5))

# Get the magma colormap
cmap = plt.get_cmap('magma')
cmap2 = plt.get_cmap('Blues')

initial_rad = "same_C" # "same_n", "same_C", "log_space_same_C", "
epilson = "0.0001" # 0.01, 0.001, 0.0001, 1.0e-5, 1.0e-8
aspect_ratio = 50

if float(epilson) < 1.0e-5:
    time_scale_factor = (3/(2*aspect_ratio))**(2/3)
else:
    time_scale_factor = (3/(2*aspect_ratio))**(1)

axlabels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)']

combo = 2

if combo == 0: # plot p = 1, c = 1 and p = 2, c = 1
    SD = np.load("/Users/cotton/Documents/Oceananigans/1DFrazil.jl/new_v_" + initial_rad + "_just_collisions_epsilon" + epilson + "_200_sum_nj.npy")
    SD2 = np.load("/Users/cotton/Documents/Oceananigans/1DFrazil.jl/new_v_" + initial_rad + "_just_collisions_epsilon" + epilson + "_200.npy")
elif combo == 1: # plot p = 1, c = 1 and p = 2, c = 3
    SD = np.load("/Users/cotton/Documents/Oceananigans/1DFrazil.jl/new_v_" + initial_rad + "_just_collisions_epsilon" + epilson + "_200_new_spherical_sum_nj.npy")
    SD2 = np.load("/Users/cotton/Documents/Oceananigans/1DFrazil.jl/new_v_" + initial_rad + "_just_collisions_epsilon" + epilson + "_200.npy")
elif combo == 2: # plot p = 2, c = 3 and p = 2, c = 3, r = 2
    SD2 = np.load("/Users/cotton/Documents/Oceananigans/1DFrazil.jl/new_v_" + initial_rad + "_just_collisions_epsilon" + epilson + "_200_new_spherical_sum_nj_r_max_redistribute.npy")
    SD = np.load("/Users/cotton/Documents/Oceananigans/1DFrazil.jl/new_v_" + initial_rad + "_just_collisions_epsilon" + epilson + "_200_new_spherical_sum_nj_redistribute.npy")

numSizeClasses = 200

if initial_rad == "log_space_same_C":
    Rs = np.logspace(np.log10(0.01*1e-3), np.log10(2*1e-3), 200)
else:
    Rs = np.linspace(0.01, 2, numSizeClasses) * 1e-3


aspect_ratio = 50
cfact = 2*np.pi*Rs**3/(aspect_ratio)

SD_to_plot = [1, 2, 5, 10, 50, 100, 200]
n_colors = len(SD_to_plot) + 1
colours = [cmap(i / (n_colors - 1)) for i in range(n_colors)]

n1max = 1.298846543545204 * 1e13
fractions = np.array([0.1, 0.9, 0.99, 0.999])
n_colors2 = len(fractions) + 1
colours2 = [cmap2(i / (n_colors2 - 1)) for i in range(n_colors2)]

if combo == 0:
    plot_data(axm, SD, 1, "CD", 1)
    plot_data(axm, SD2, 2, "CP", 1)
elif combo == 1:
    plot_data(axm, SD, 1, "SD", 1)
    plot_data(axm, SD2, 2, "CP", 1)
elif combo == 2:
    plot_data(axm, SD, 1, "SDR", time_scale_factor)
    plot_data(axm, SD2, 2, "SD", 1)

tend = np.max([time_scale_factor*SD[0, np.argmin(np.abs(SD[1, :] - (1-1e-6)*n1max))], SD2[0, np.argmin(np.abs(SD2[1, :] - (1-1e-6)*n1max))]])

axm["n"].set_xlim([0, tend])    
axm["c"].set_xlim([0, tend])
axm["t"].set_xlim([0, tend])  
axm["p"].set_xlim([Rs[1]*1000, Rs[-1]*1000])  


ymaxp = max(axm["p"].get_ylim()) * 1.05
yminp = min(axm["p"].get_ylim())
axm["p"].set_ylim([yminp, ymaxp])        


axm["n"].set_ylim([1, 2e13])    
axm["c"].set_ylim([1e-10, 0.002])  
axm["t"].set_ylim([-1, 1]) 
axm["t"].set_yticks([])    
axm["t"].set_yticklabels([])    
axm["n"].legend(frameon = True, loc = "upper right", bbox_to_anchor = (1, 0.99))
axm["c"].legend(frameon = True, loc = "upper right", bbox_to_anchor = (1, 0.99))
axm["p"].legend(frameon = True)
#axm["t"].legend(frameon = True)

axm["c"].set_xlabel("time (s)")
axm["n"].set_xlabel("time (s)")
axm["n"].set_ylabel(r"$n_i$")
axm["c"].set_ylabel(r"$C_i$")

if initial_rad == "same_n":
    axm["p"].set_ylabel(r"$n(r)$")
else:
    axm["p"].set_ylabel(r"$C(r)$")
    
axm["p"].set_xlabel(r"r (mm)")
axm["t"].set_xlabel("time (s)")

axm["n"].text(0.01, 0.95, axlabels[0],
    transform=axm["n"].transAxes,
    fontsize=12,
    verticalalignment='bottom',
    horizontalalignment='left')
axm["c"].text(0.01, 0.95, axlabels[1],
    transform=axm["c"].transAxes,
    fontsize=12,
    verticalalignment='bottom',
    horizontalalignment='left')
axm["p"].text(0.01, 0.93, axlabels[2],
    transform=axm["p"].transAxes,
    fontsize=12,
    verticalalignment='bottom',
    horizontalalignment='left')
axm["t"].text(0.01, 0.65, axlabels[3],
    transform=axm["t"].transAxes,
    fontsize=12,
    verticalalignment='bottom',
    horizontalalignment='left')

SD_times, SD2_times = find_time_scale_factor(SD, SD2)
print(SD2_times/SD_times *1/time_scale_factor)

test  = SD2_times/SD_times

# Collect legend handles and labels from one or more axes
handles, labels = [], []
h, l = axm["t"].get_legend_handles_labels()
handles.extend(h)
labels.extend(l)

# Add shared legend below all axes
fig.legend(handles, labels,
           loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.02), frameon = False)

# Adjust layout to make room for the legend
fig.tight_layout(rect=[0, 0.02, 1, 1])  # Reserve space at bottom

if combo == 0:
    plt.savefig("21_" + initial_rad + "_eps_" + epilson + ".pdf")
elif combo == 1:
    plt.savefig("23_" + initial_rad + "_eps_" + epilson + ".pdf")
elif combo == 2:
    plt.savefig("232_" + initial_rad + "_eps_" + epilson + "tscale.pdf")
