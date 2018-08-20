# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 22:30:43 2016
@author: Brendan Gallagher

"""


import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import f as ff
from matplotlib import cm
import yaml
import sunpy.cm
import sys
import os 

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27

with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

directory = cfg['processed_dir']
date = cfg['date']
wavelength = cfg['wavelength']
savefig = cfg['save_fig']

# 11-param-list
titles = [r'Power Law Slope-Coefficient [flux] - A', r'Power Law Index n', 
          r'Power Law Tail - C', r'Lorentzian Amplitude [flux] - α', 
          r'Lorentz. Loc. β [min]', r'Lorentzian Width - σ', 'F-Statistic', 
          r'Lorentzian Amplitude Scaled - α', r'$r$-Value: Correlation Coefficient', 
          r'Rollover Period $T_r$ [min]', r'$\chi^2$']
names = ['slope_coeff', 'index', 'tail', 'lorentz_amp', 'lorentz_loc', 
         'lorentz_wid', 'f_test', 'lorentz_amp_scaled', 'r_value', 'roll_freq', 'chisqr']

# load parameter/heatmap array 
h_map = np.load('%s/param.npy' % directory)  

# generate p-value heatmap + masked Lorentzian component heatmaps
dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
p_val = ff.sf(h_map[6], dof1, dof2)

mask_thresh = 0.005  # significance threshold
   
p_mask = np.copy(p_val)
amp_mask = np.copy(h_map[3])
loc_mask = np.copy(h_map[4])
wid_mask = np.copy(h_map[5])    

# mask the Lorenztian component arrays with NaNs if above threshold 
p_mask[p_val > mask_thresh] = np.NaN
amp_mask[p_val > mask_thresh] = np.NaN
loc_mask[p_val > mask_thresh] = np.NaN
wid_mask[p_val > mask_thresh] = np.NaN    

# determine percentage of region masked 
count = np.count_nonzero(np.isnan(p_mask))   
total_pix = p_val.shape[0]*p_val.shape[1]
mask_percent = ((np.float(count))/total_pix)*100

# convert Lorentzian location to minutes
h_map[4] = (1./np.exp(h_map[4]))/60.               
loc_mask = (1./np.exp(loc_mask))/60.  

plots = [p_mask, amp_mask, loc_mask, wid_mask]

if h_map.shape[2] > h_map.shape[1]:
    aspect_ratio = float(h_map.shape[2]) / float(h_map.shape[1])
    fig_height = 10
    fig_width = 10*aspect_ratio
    
else:
    aspect_ratio = float(h_map.shape[1]) / float(h_map.shape[2])
    fig_width = 10
    fig_height = 10*aspect_ratio


for i in range(h_map.shape[0]):
    
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca() 
    plt.title(r'%s' % (titles[i]), y = 1.02, fontsize=font_size)
    
    param = h_map[i]
    
    # specify discrete colorscale with 10 intervals
    if i == 4:
        cmap = cm.get_cmap('jet_r', 10)
    else: 
        cmap = cm.get_cmap('jet', 10)   
     
    pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
    pNaN = pflat[~np.isnan(pflat)]
    h_min = np.percentile(pNaN,1)
    h_max = np.percentile(pNaN,99)

    # specify colorbar ticks to be at boundaries of segments
    h_range = np.abs(h_max-h_min)
    h_step = h_range / 10.
    c_ticks = np.zeros((11))
    for h in range(11):
        c_ticks[h] = h_min + h_step*h 
        
    im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=h_min, vmax=h_max)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    ax.tick_params(axis='both', which='major', pad=10)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax, format='%0.2f')
    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    cbar.set_ticks(c_ticks)
    
    if savefig == True:
        plt.savefig('%s/Figures/%s.pdf' % (directory, names[i]), 
                    format='pdf', bbox_inches='tight')
    
    
    if i == 3 or i == 4 or i == 5:   
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca() 
        
        plt.title(r'%s; p < %0.3f | f$_{masked}$ = %0.1f%s' % 
                  (titles[i], mask_thresh, mask_percent, '%'), 
                  y = 1.02, fontsize=font_size)

        if i == 4:
            cmap = cm.get_cmap('jet_r', 10)
        else:
            cmap = cm.get_cmap('jet', 10)                    
            
        im = ax.imshow(np.flipud(plots[i-2]), cmap=cmap, vmin=h_min, vmax=h_max)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
        ax.tick_params(axis='both', which='major', pad=10)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax, format='%0.2f')
        cbar.ax.tick_params(labelsize=font_size, pad=5) 
        cbar.set_ticks(c_ticks)

        if savefig == True:
            plt.savefig('%s/Figures/%s_mask_%i.pdf' % 
                        (directory, names[i], (1./mask_thresh)), 
                        format='pdf', bbox_inches='tight')


# generate visual images
vis = np.load('%s/visual.npy' % directory)  
vis = vis[1:-1,1:-1]  # make same size as heatmaps (if using 3x3 averaging)
          
v_min = np.percentile(vis,1)
v_max = np.percentile(vis,99) 

fig = plt.figure(figsize=(fig_width,fig_height))
ax = plt.gca()
plt.title('Visual Average', y = 1.02, fontsize=font_size)
im = ax.imshow(np.flipud(vis), cmap='sdoaia%i' % wavelength, vmin=v_min, vmax=v_max)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
cbar.ax.tick_params(labelsize=font_size, pad=5) 

if savefig == True:
    plt.savefig('%s/Figures/visual_average.pdf' % directory, 
                format='pdf', bbox_inches='tight')
    
scriptName = os.path.splitext(os.path.basename(sys.argv[0]))[0]
  
#with open('%s_%i_region_details.txt' % (date, wavelength), 'w') as file:
with open('log.txt', 'a+') as file:
    file.write("%s: Generate Parameter Heatmaps" % scriptName + "\n")
    file.write("--------------------------------------" + "\n")
    file.write("Parameter heatmap figures saved: %s" % savefig + "\n\n")
    file.write("========================================" + "\n")
    file.write("++++++++++++++++++++++++++++++++++++++++" + "\n")
    file.write("========================================" + "\n\n")