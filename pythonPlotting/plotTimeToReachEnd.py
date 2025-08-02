#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 18:04:43 2025

@author: cotton

Plot time to reach end

"""

# get the time to reach x %

import matplotlib.pyplot as plt
import numpy as np

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

axlabels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)', '(p)']

fig, axm = plt.subplot_mosaic(mosaic= """
nN
cC
lL
""", sharex = True, layout="tight", figsize=(8, 10))



def get_times(epsilons, string_end, initial_rad, n1val):
    times = np.zeros(len(epsilons))
    rbars = np.zeros(len(epsilons))
    for indx, epsilon in enumerate(epsilons):
        ns = np.load("/Users/cotton/Documents/Oceananigans/1DFrazil.jl/new_v_" + initial_rad + "_just_collisions_epsilon" + epsilon + string_end)
        tindx = np.argmin(np.abs(ns[1, :]- n1val))
        times[indx] = ns[0, tindx]
        rbars[indx] = np.mean(ns[2:, tindx]*Rs[1:]/sum(ns[2:, tindx]))*1000

    
    return times, rbars

def plot_times_firstChapter():
    for i, n1val in enumerate(fractions*n1max):
        frac = fractions[i]
        label_str = r"$\alpha = " + str(frac) + "$"

    
    
        initial_rad = "same_n" # "log_space_same_C", "
        # time for CP
        string_end = "_200.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["n"].loglog(epsilonsn, times, color = colours2[i+1], ls = "--", label = label_str)
        axm["N"].semilogx(epsilonsn, rbars, color = colours2[i+1], ls = "--")
    
        # time for SD
        string_end = "_200_new_spherical_sum_nj.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["n"].loglog(epsilonsn, times, color = colours2[i+1])
        axm["N"].semilogx(epsilonsn, rbars, color = colours2[i+1])
        
        
        
        initial_rad = "same_C" # "log_space_same_C", "
        # time for CP
        string_end = "_200.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["c"].loglog(epsilonsn, times, color = colours2[i+1], ls = "--", label = "CP")
        axm["C"].semilogx(epsilonsn, rbars, color = colours2[i+1], ls = "--")
    
        # time for SD
        string_end = "_200_new_spherical_sum_nj.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["c"].loglog(epsilonsn, times, color = colours2[i+1], label = "SD")
        axm["C"].semilogx(epsilonsn, rbars, color = colours2[i+1])
        
        initial_rad = "log_space_same_C" # "log_space_same_C", "
        # time for CP
        string_end = "_200.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["l"].loglog(epsilonsn, times, color = colours2[i+1], ls = "--", label = "CP")
        axm["L"].semilogx(epsilonsn, rbars, color = colours2[i+1], ls = "--")
    
        # time for SD
        string_end = "_200_new_spherical_sum_nj.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["l"].loglog(epsilonsn, times, color = colours2[i+1], label = "SD")
        axm["L"].semilogx(epsilonsn, rbars, color = colours2[i+1])
        
    axm["n"].loglog(10, times[0], color = "k", ls = "-", label = "SD")
    axm["n"].loglog(10, times[0], color = "k", ls = "--", label = "CP")
        
    axm["n"].legend()

def plot_times_secondChapter():
    for i, n1val in enumerate(fractions*n1max):
        frac = fractions[i]
        label_str = r"$\alpha = " + str(frac) + "$"

    
    
        initial_rad = "same_n" # "log_space_same_C", "
        # time for CP
        string_end = "_200_new_spherical_sum_nj.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["n"].loglog(epsilonsn, times, color = colours2[i+1], ls = "--", label = label_str)
        axm["N"].semilogx(epsilonsn, rbars, color = colours2[i+1], ls = "--")
    
        # time for SD
        string_end = "_200_new_spherical_sum_nj_redistribute.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["n"].loglog(epsilonsn, times, color = colours2[i+1])
        axm["N"].semilogx(epsilonsn, rbars, color = colours2[i+1])
        
        
        
        initial_rad = "same_C" # "log_space_same_C", "
        # time for CP
        string_end = "_200_new_spherical_sum_nj.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["c"].loglog(epsilonsn, times, color = colours2[i+1], ls = "--", label = "CP")
        axm["C"].semilogx(epsilonsn, rbars, color = colours2[i+1], ls = "--")
    
        # time for SD
        string_end = "_200_new_spherical_sum_nj_redistribute.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["c"].loglog(epsilonsn, times, color = colours2[i+1], label = "SD")
        axm["C"].semilogx(epsilonsn, rbars, color = colours2[i+1])
        
        initial_rad = "log_space_same_C" # "log_space_same_C", "
        # time for CP
        string_end = "_200_new_spherical_sum_nj.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["l"].loglog(epsilonsn, times, color = colours2[i+1], ls = "--", label = "CP")
        axm["L"].semilogx(epsilonsn, rbars, color = colours2[i+1], ls = "--")
    
        # time for SD
        string_end = "_200_new_spherical_sum_nj_redistribute.npy"
        times, rbars = get_times(epsilons, string_end, initial_rad, n1val)
        axm["l"].loglog(epsilonsn, times, color = colours2[i+1], label = "SD")
        axm["L"].semilogx(epsilonsn, rbars, color = colours2[i+1])
        
    axm["n"].loglog(10, times[0], color = "k", ls = "-", label = "SDR")
    axm["n"].loglog(10, times[0], color = "k", ls = "--", label = "SD")
        
    axm["n"].legend()

def add_lettering():
    axm["n"].text(0.01, 0.92, axlabels[0],
        transform=axm["n"].transAxes,
        fontsize=12,
        verticalalignment='bottom',
        horizontalalignment='left')
    axm["N"].text(0.01, 0.92, axlabels[1],
        transform=axm["N"].transAxes,
        fontsize=12,
        verticalalignment='bottom',
        horizontalalignment='left')
    axm["c"].text(0.01, 0.92, axlabels[2],
        transform=axm["c"].transAxes,
        fontsize=12,
        verticalalignment='bottom',
        horizontalalignment='left')
    axm["C"].text(0.01, 0.92, axlabels[3],
        transform=axm["C"].transAxes,
        fontsize=12,
        verticalalignment='bottom',
        horizontalalignment='left')
    axm["l"].text(0.01, 0.92, axlabels[4],
        transform=axm["l"].transAxes,
        fontsize=12,
        verticalalignment='bottom',
        horizontalalignment='left')
    axm["L"].text(0.01, 0.92, axlabels[5],
        transform=axm["L"].transAxes,
        fontsize=12,
        verticalalignment='bottom',
        horizontalalignment='left')

def expand_limits(keys):
    for key in keys:
        ymax = max(axm[key].get_ylim()) * 1.3
        ymin = min(axm[key].get_ylim())
        axm[key].set_ylim([ymin, ymax])

cmap2 = plt.get_cmap('Blues')
fractions = np.array([0.1, 0.9, 0.99, 0.999])
n_colors2 = len(fractions) + 1
colours2 = [cmap2(i / (n_colors2 - 1)) for i in range(n_colors2)]    
n1max = 1.298846543545204 * 1e13

epsilons = np.array(["0.01", "0.001", "0.0001", "1.0e-5", "1.0e-8"])
epsilonsn = epsilons.astype(float)
numSizeClasses = 200
Rs = np.linspace(0.01, 2, numSizeClasses) * 1e-3

#plot_times_firstChapter()
plot_times_secondChapter()

axm["n"].set_title("same number, linearly spaced")
axm["c"].set_title("same concentration, linearly spaced")
axm["l"].set_title("same concentration, logarithmically spaced")
axm["N"].set_title("same number, linearly spaced")
axm["C"].set_title("same concentration, linearly spaced")
axm["L"].set_title("same concentration, logarithmically spaced")

axm["n"].set_xlim([np.min(epsilonsn), np.max(epsilonsn)]) 
axm["n"].set_ylabel(r"$t \,(n_1 = \alpha n_{max})$")
axm["N"].set_ylabel(r"$\bar{r} \, (n_1 = \alpha n_{max})$")
axm["c"].set_ylabel(r"$t \,(n_1 = \alpha n_{max})$")
axm["C"].set_ylabel(r"$\bar{r} \, (n_1 = \alpha n_{max})$")
axm["l"].set_ylabel(r"$t \,(n_1 = \alpha n_{max})$")
axm["L"].set_ylabel(r"$\bar{r} \, (n_1 = \alpha n_{max})$")
axm["l"].set_xlabel(r"$\epsilon (m^2 s^{-3})$")
axm["L"].set_xlabel(r"$\epsilon (m^2 s^{-3})$")

add_lettering()
expand_limits(["n", "c", "l"])

# time for SDR
#string_end = "_200_new_spherical_sum_nj_redistribute.npy"
#axm["n"].semilogx(epsilonsn, times, color = CB_color_cycle[2], ls = ":")
#axm["N"].semilogx(epsilonsn, rbars, color = CB_color_cycle[2], ls = ":")