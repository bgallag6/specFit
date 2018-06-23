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
import sunpy

   
# define Gaussian-fitting function
def Gauss(f, P, fp, fw):
    return P*np.exp(-0.5*((f-fp)/fw)**2)

with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

directory = cfg['temp_dir']
date = cfg['date']
wavelength = cfg['wavelength']
savefig = False

  
# 11-param-list
titles = [r'Power Law Slope-Coefficient [flux] - A', r'(b) Power Law Index n', r'Power Law Tail - C', r'Gaussian Amplitude [flux] - α', r'(c) Gauss. Loc. β [min]', r'Gaussian Width - σ', 'F-Statistic', r'Gaussian Amplitude Scaled - α', r'$r$-Value: Correlation Coefficient', r'(d) Rollover Period $T_r$ [min]', r'$\chi^2$']
names = ['slope_coeff', 'index', 'tail', 'lorentz_amp', 'lorentz_loc', 'lorentz_wid', 'f_test', 'lorentz_amp_scaled', 'r_value', 'roll_freq', 'chisqr']

# load parameter array and visual images from file tree structure 
heatmaps = np.load('%s/Processed/Output/%s/%i/param.npy' % (directory, date, wavelength))
visual = np.load('%s/Processed/Output/%s/%i/visual.npy'% (directory, date, wavelength))  

visual = visual[1:-1,1:-1]  # to make same size as heatmaps (if using 3x3 pixel box averaging)
h_map = heatmaps    

plt.rcParams["font.family"] = "Times New Roman"
font_size = 27  # set the font size to be used for all text - titles, tick marks, text, labels    

xdim = int(np.floor(h_map.shape[2]/100))
ydim = int(np.floor(h_map.shape[1]/100))

x_ticks = [100*i for i in range(xdim+1)]
y_ticks = [100*i for i in range(ydim+1)]

x_ind = x_ticks
y_ind = y_ticks    

# generate p-value heatmap + masked Gaussian component heatmaps
df1, df2 = 3, 6  # degrees of freedom for model M1, M2
p_val = ff.sf(h_map[6], df1, df2)

mask_thresh = 0.005  # significance threshold - masked above this value
   
p_mask = np.copy(p_val)
amp_mask = np.copy(h_map[3])
loc_mask = np.copy(h_map[4])
wid_mask = np.copy(h_map[5])    

# mask the Gaussian component arrays with NaNs if above threshold 
p_mask[p_val > mask_thresh] = np.NaN  # for every element in p_mask, if the corresponding element in p_val is greater than the threshold, set that value to NaN
amp_mask[p_val > mask_thresh] = np.NaN
loc_mask[p_val > mask_thresh] = np.NaN
wid_mask[p_val > mask_thresh] = np.NaN    

# determine percentage of region masked 
count = np.count_nonzero(np.isnan(p_mask))   
total_pix = p_val.shape[0]*p_val.shape[1]
mask_percent = ((np.float(count))/total_pix)*100

h_map[4] = (1./np.exp(h_map[4]))/60.               
loc_mask = (1./np.exp(loc_mask))/60.  # convert Gaussian location to minutes
plots = [p_mask, amp_mask, loc_mask, wid_mask]  # make array of masked plots to iterate over

if h_map.shape[2] > h_map.shape[1]:
    aspect_ratio = float(h_map.shape[2]) / float(h_map.shape[1])
    fig_height = 10
    fig_width = 10*aspect_ratio
    
else:
    aspect_ratio = float(h_map.shape[1]) / float(h_map.shape[2])
    fig_width = 10
    fig_height = 10*aspect_ratio  # works better for 20130626


for i in range(h_map.shape[0]):
    
    fig = plt.figure(figsize=(fig_width,fig_height))
    ax = plt.gca()  # get current axis -- to set colorbar 
    plt.title(r'%s' % (titles[i]), y = 1.02, fontsize=font_size)
    
    param = h_map[i]

    if i == 4:
        cmap = cm.get_cmap('jet_r', 10)
    else: 
        cmap = cm.get_cmap('jet', 10)  # specify discrete colorscale with 10 intervals 
     
    pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
    pNaN = pflat[~np.isnan(pflat)]
    h_min = np.percentile(pNaN,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max = np.percentile(pNaN,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)

    # specify colorbar ticks to be at boundaries of segments
    h_range = np.abs(h_max-h_min)
    h_step = h_range / 10.
    c_ticks = np.zeros((11))
    for h in range(11):
        c_ticks[h] = h_min + h_step*h 
        
    im = ax.imshow(np.flipud(h_map[i]), cmap = cmap, vmin=h_min, vmax=h_max)
    plt.xticks(x_ticks,x_ind,fontsize=font_size)
    plt.yticks(y_ticks,y_ind,fontsize=font_size)
    ax.tick_params(axis='both', which='major', pad=10)
    divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax, format='%0.2f')

    cbar.ax.tick_params(labelsize=font_size, pad=5) 
    cbar.set_ticks(c_ticks)
    
    if savefig == True:
        plt.savefig('%s/Processed/Output/%s/%i/Figures/%s_%i_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf', bbox_inches='tight')
    
    
    if i == 3 or i == 4 or i == 5:   
        fig = plt.figure(figsize=(fig_width,fig_height))
        ax = plt.gca()  # get current axis -- to set colorbar 
        
        plt.title(r'%s; p < %0.3f | f$_{masked}$ = %0.1f%s' % (titles[i], mask_thresh, mask_percent, '%'), y = 1.02, fontsize=font_size)

        if i == 4:
            cmap = cm.get_cmap('jet_r', 10)
        else:
            cmap = cm.get_cmap('jet', 10)                    
            
        im = ax.imshow(np.flipud(plots[i-2]), cmap = cmap, vmin=h_min, vmax=h_max)
        plt.xticks(x_ticks, x_ind, fontsize=font_size)
        plt.yticks(y_ticks, y_ind, fontsize=font_size)
        ax.tick_params(axis='both', which='major', pad=10)
        divider = make_axes_locatable(ax)  # set colorbar to heatmap axis
        cax = divider.append_axes("right", size="3%", pad=0.07)
        cbar = plt.colorbar(im,cax=cax, format='%0.2f')
        cbar.ax.tick_params(labelsize=font_size, pad=5) 
        cbar.set_ticks(c_ticks)

        if savefig == True:
            plt.savefig('%s/Processed/Output/%s/%i/Figures/%s_%i_%s_mask_%i.pdf' % (directory, date, wavelength, date, wavelength, names[i], (1./mask_thresh)), format='pdf', bbox_inches='tight')
        
    
    flat_param = np.reshape(h_map[i], (h_map[i].shape[0]*h_map[i].shape[1]))
    
    # calculate some statistics
    sigma = np.std(flat_param)   
    
    fig = plt.figure(figsize=(fig_width+1,fig_height))
    plt.title('%s' % (titles[i]), y = 1.02, fontsize=font_size)  # no date / wavelength
    plt.xlabel('%s' % titles[i], fontsize=font_size, labelpad=10)
    plt.ylabel('Bin Count', fontsize=font_size, labelpad=10)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    plt.xlim(h_min, h_max)
    y, x, _ = plt.hist(flat_param, bins=200, range=(h_min, h_max))
    
    #n, bins, patches = plt.hist(flat_param, bins=200, range=(h_min, h_max))
    n=y[1:-2]
    bins=x[1:-2]
    elem = np.argmax(n)
    bin_max = bins[elem]
    plt.ylim(0, y.max()*1.1)      
    
    
    plt.vlines(bin_max, 0, y.max()*1.1, color='black', linestyle='dotted', linewidth=2., label='mode=%0.4f' % bin_max)  
    plt.vlines(0, 0, y.max()*1.1, color='white', linestyle='dashed', linewidth=1.5, label='sigma=%0.4f' % sigma)
    legend = plt.legend(loc='upper right', prop={'size':20}, labelspacing=0.35)
    for label in legend.get_lines():
        label.set_linewidth(2.0)  # the legend line width
    
    if savefig == True:
        plt.savefig('%s/Processed/Output/%s/%i/Figures/%s_%i_Histogram_%s.pdf' % (directory, date, wavelength, date, wavelength, names[i]), format='pdf', bbox_inches='tight')


# generate visual images
vis = visual
          
v_min = np.percentile(vis,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
v_max = np.percentile(vis,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)  

fig = plt.figure(figsize=(fig_width,fig_height))
ax = plt.gca()
#ax = plt.subplot2grid((1,31),(0, 0), colspan=30, rowspan=1)  #to substitute for colorbar space
#plt.subplots_adjust(right=0.875)  #to substitute for colorbar space
plt.title('Visual Average', y = 1.02, fontsize=font_size, fontname="Times New Roman")  # no date / wavelength
#plt.title('(f) %i $\AA$' % wavelength, y = 1.02, fontsize=font_size, fontname="Times New Roman")  # no date / wavelength
#im = ax.imshow(h_map[i], vmin=vmin[i], vmax=vmax[i])
#im = ax.imshow(vis[i], cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)   
im = ax.imshow(np.flipud(vis), cmap='sdoaia%i' % wavelength, vmin = v_min, vmax = v_max)
plt.xticks(x_ticks,x_ind,fontsize=font_size, fontname="Times New Roman")
plt.yticks(y_ticks,y_ind,fontsize=font_size, fontname="Times New Roman")
ax.tick_params(axis='both', which='major', pad=10)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.07)
cbar = plt.colorbar(im,cax=cax)
cbar.ax.tick_params(labelsize=font_size, pad=5) 

if savefig == True:
    plt.savefig('%s/DATA/Output/%s/%i/Figures/%s_%i_visual_average.pdf' % (directory, date, wavelength, date, wavelength), format='pdf', bbox_inches='tight')