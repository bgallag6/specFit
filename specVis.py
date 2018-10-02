# -*- coding: utf-8 -*-
"""
Created on Thu May 24 08:25:54 2018

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit as Fit
from pylab import axvline
import sunpy
import sunpy.cm
from scipy import fftpack
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Button, Slider
from scipy.stats import f as ff
from scipy.stats.stats import pearsonr
import os
import yaml
from specModel import M1, M2, m2

    
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
def update(val):
    global mask_val
    mask_val = slid_mask.val
    return mask_val

def colorBar():
    global cax1
    global cbar1
    # design colorbar for heatmaps
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="3%", pad=0.07)
    cbar1 = plt.colorbar(im,cax=cax1)
    cbar1.ax.tick_params(labelsize=13, pad=3)   
    plt.colorbar(im,cax=cax1)
    plt.draw()

def plotMap(p):
    global cbar1
    global im
    cbar1.remove()
    im.remove()
    param = h_map[p]
    pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
    pNaN = pflat[~np.isnan(pflat)]
    h_min = np.percentile(pNaN,1)
    h_max = np.percentile(pNaN,99)
    if p == 4:
        c_map = 'jet_r'
    else:
        c_map = 'jet'
    im = ax1.imshow(param, cmap=c_map, interpolation='nearest', vmin=h_min, 
                    vmax=h_max, picker=True)
    ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[p]), 
                  y = 1.01, fontsize=17)
    colorBar()
    
def plotMask(p):
    global mask_val
    global cbar1
    global im
    cbar1.remove()
    im.remove()  
    
    param = h_map[p]
    pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
    pNaN = pflat[~np.isnan(pflat)]
    h_min = np.percentile(pNaN,1)
    h_max = np.percentile(pNaN,99)
    
    # generate p-value heatmap + masked Lorentzian component heatmaps
    dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
    p_val = ff.sf(h_map[6], dof1, dof2)
    param_mask = np.copy(param) 
    param_mask[p_val > mask_val] = np.NaN
    
    # determine percentage of region masked 
    count = np.count_nonzero(np.isnan(param_mask))   
    total_pix = p_val.shape[0]*p_val.shape[1]
    mask_percent = ((np.float(count))/total_pix)*100
    
    if p == 4:
        c_map = 'jet_r'
    else:
        c_map = 'jet'
    im = ax1.imshow(param_mask, cmap=c_map, interpolation='nearest', 
                    vmin=h_min, vmax=h_max, picker=True)
    ax1.set_title(r'%s: %i $\AA$ | %s | $f_{masked}$ = %0.1f%s' % 
                  (date_title, wavelength, titles[p], mask_percent, '%'), 
                  y=1.01, fontsize=17)
    colorBar()
        
def histMask(p):
    global mask_val
    param = h_map[p]
    
    # generate p-value heatmap + masked Lorentzian component heatmaps
    dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
    p_val = ff.sf(h_map[6], dof1, dof2)
    param_mask = np.copy(param) 
    param_mask[p_val > mask_val] = 0.
    param1d = np.reshape(param_mask, (param_mask.shape[0]*param_mask.shape[1]))
    pmask = param1d[param1d != 0]
    pmask = pmask[~np.isnan(pmask)]
    ax2.set_title('Histogram: %s | Masked' % titles[marker], y = 1.01, fontsize=17)
    ax2.hist(pmask, bins=25, edgecolor='black')

    
class Index(object):
    ind = 0
         
    def coeff(self, event):
        global marker
        marker = 0
        plotMap(marker)  # could just do this
        return marker       

    def index(self, event):
        global marker
        marker = 1
        plotMap(marker)
        return marker     
        
    def tail(self, event):  # meh, should probably fix this
        global marker
        marker = 2
        plotMap(marker)
        return marker
        
    def lorentz_amp(self, event):
        global marker
        marker = 3
        plotMap(marker)
        return marker
    
    def lorentz_loc(self, event):
        global marker
        marker = 4
        plotMap(marker)
        return marker
        
    def lorentz_wid(self, event):
        global marker
        marker = 5
        plotMap(marker)
        return marker
        
    def fstat(self, event):
        global marker
        marker = 6
        plotMap(marker)
        return marker
        
    def visual(self, event):
        global cbar1
        global im
        cbar1.remove()
        im.remove()
        param = vis
        h_min = np.percentile(param,1)
        h_max = np.percentile(param,99)
        im = ax1.imshow(param, cmap='sdoaia%i' % wavelength, interpolation='nearest', 
                        vmin=h_min, vmax=h_max, picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[7]), 
                      y = 1.01, fontsize=17)
        colorBar()
    
    def hist(self, event):
        global marker
        global toggle2
        ax2.clear()
        if toggle2 == 0:
            ax2.set_title('Histogram: %s' % titles[marker], y = 1.01, fontsize=17)
            pflat = np.reshape(h_map[marker], (h_map[marker].shape[0]*h_map[marker].shape[1]))
            pNaN = pflat[~np.isnan(pflat)]
            ax2.hist(pNaN, bins=25, edgecolor='black')
        elif toggle2 == 1:
            histMask(marker)
        plt.draw()

        
    def mask(self, event):
        global toggle2
        global marker
        if toggle2 == 0:
            toggle2 = 1
            plotMask(marker)
        elif toggle2 == 1:
            toggle2 = 0
            plotMap(marker)              
        return toggle2
        
    def saveFig(self, event):
        global count

        outdir = 'C:/Users/Brendan/Desktop/Tool_Figures'

        if not os.path.exists(os.path.dirname(outdir)):
            try:
                print("Specified directory not found.")
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST: raise
        else:
            plt.savefig('%s/%s_%i_%i.pdf' % (outdir, date, wavelength, count), 
                        bbox_inches='tight')
        count += 1
        return count
  

# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy
    ixx, iyy = event.xdata, event.ydata
    if ixx > 1. and iyy > 1.:
        ax2.clear()
        del ax1.collections[:]
        plt.draw()
        print ('x = %i, y = %i' % ( ixx, iyy))  # print location of pixel
        ix = int(ixx)
        iy = int(iyy)
        
        s = np.zeros((spectra.shape[2]))
    
        s[:] = spectra[iy][ix][:]
        
        # use 3x3 pixel-box std.dev. or adhoc method for fitting uncertainties
        if spec_unc == 'stddev':
            ds = np.zeros((spectra.shape[2]))
            ds[:] = stdDev[iy][ix][:]
        elif spec_unc == 'adhoc':
            ds = ds0
     
        ## fit data to combined power law plus gaussian component model
        
        try:
            nlfit_l, nlpcov_l = Fit(M1, freqs, s, p0=M1_guess, 
                                    bounds=(M1_low, M1_high), sigma=ds, 
                                    method='dogbox')
          
        except RuntimeError: print("Error M1 - curve_fit failed")
        except ValueError: print("Error M1 - inf/NaN")
        
        A, n, C = nlfit_l  # unpack fitting parameters
         
       
        try:                                                          
            nlfit_gp, nlpcov_gp = Fit(M2, freqs, s, p0=M2_guess, 
                                      bounds=(M2_low,M2_high), sigma=ds, 
                                      method='dogbox', max_nfev=3000)
        
        except RuntimeError: print("Error M2 - curve_fit failed")
        except ValueError: print("Error M2 - inf/NaN")

        
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
        #print nlfit_gp
               
        try:           
            nlfit_gp2, nlpcov_gp2 = Fit(M2, freqs, s, p0=nlfit_gp, 
                                        bounds=(M2_low, M2_high), sigma=ds, 
                                        max_nfev=3000)    
           
        except RuntimeError: print("Error M2 - curve_fit failed")
        except ValueError: print("Error M2 - inf/NaN")
        
        A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack model parameters     
        #print nlfit_gp2
                       
        # create models from parameters    
        m1_fit = M1(freqs, A, n, C)    
        lorentz = m2(freqs,P22,fp22,fw22)
        m2_fit2 = M2(freqs, A22,n22,C22,P22,fp22,fw22) 
        m1_fit2 = M1(freqs, A22,n22,C22)      
        
        residsM1 = (s - m1_fit)
        chisqrM1 =  ((residsM1/ds)**2).sum() 
        
        residsM22 = (s - m2_fit2)
        chisqrM22 = ((residsM22/ds)**2).sum()
        redchisqrM22 = ((residsM22/ds)**2).sum()/float(freqs.size-6) 
              
        f_test2 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(freqs.size-6))
        
        dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
        p_val = ff.sf(f_test2, dof1, dof2)
        
        # extract the lorentzian amplitude scaling factor
        amp_scale2 = M1(np.exp(fp22), A22, n22, C22)  

        rval = pearsonr(m2_fit2, s)[0]
        
        rollover = (1. / ((C22 / A22)**(-1. / n22))) / 60.
        
        fwhm = (1. / (np.exp(fp22+fw22) - np.exp(fp22-fw22))) / 60.
        
        pLoc = (1./np.exp(fp22))/60.
        
        ax2.loglog(freqs, s, 'k', linewidth=1.5)
        
        #ax2.loglog(f, m1_fit, 'r', linewidth=1.3, label='M1')
        ax2.loglog(freqs, m2_fit2, c='purple', lw=1.5, label='M2')
        ax2.loglog(freqs, m1_fit2, c='g', lw=1.5, label='M2: Power Law')
        ax2.loglog(freqs, lorentz, c='g', ls='--', lw=1.5, label='M2: Lorentz')
        
        ax2.set_xlabel('Frequency [Hz]', fontsize=fontSz-3, labelpad=5)
        ax2.set_ylabel('Power', fontsize=fontSz-3, labelpad=5)
        axvline(x=0.00333,color='k',ls='dashed', label='5 minutes')
        axvline(x=0.00555,color='k',ls='dotted', label='3 minutes')
        ax2.set_title('Spectra Fit: Pixel (%ix , %iy)' % (ix, iy), fontsize=fontSz-5)
        
        ax2.text(0.011, 10**-0.5, r'$n$ = {0:0.2f}'.format(n22), fontsize=fontSz-2)
        ax2.text(0.011, 10**-0.75, r'$\alpha$ = {0:0.2e}'.format(P22), fontsize=fontSz-2)
        ax2.text(0.011, 10**-1.00, r'$\beta$ = {0:0.1f} [min]'.format(pLoc), fontsize=fontSz-2)
        ax2.text(0.0061, 10**-1.25, r'$FWHM$ = {0:0.1f} [min]'.format(fwhm), fontsize=fontSz-2)
        legend = ax2.legend(loc='lower left', prop={'size':15}, labelspacing=0.35)
        ax2.set_xlim(10**-4.5, 10**-1.3)
        ax2.set_ylim(10**-5, 10**0)   
        ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
        for label in legend.get_lines():
                label.set_linewidth(2.0)  # the legend line width   
                
    return ix, iy

    

"""
##############################################################################
##############################################################################
"""

print("Starting specVis...", flush=True)

plt.rcParams["font.family"] = "Times New Roman"
fontSz = 20 

#directory = 'S:'
#directory = 'C:/Users/Brendan/Desktop/specFit/test'
#directory = './test/validation'
date = '20120606'
wavelength = 1600

with open('specFit_config.yaml', 'r') as stream:
    cfg = yaml.load(stream)

#directory = cfg['temp_dir']
#date = cfg['date']
#wavelength = cfg['wavelength']
mmap_spectra = cfg['mmap_spectra']
M1_low = cfg['M1_low']
M1_high = cfg['M1_high']
M2_low = cfg['M2_low']
M2_high = cfg['M2_high']
spec_unc = cfg['spec_unc']
M1_guess = cfg['M1_guess']
M2_guess = cfg['M2_guess']
directory = cfg['specVis_dir']

global spectra

if mmap_spectra == True:
    # load memory-mapped array as read-only
    spectra = np.load('%s/specCube.npy' % directory, mmap_mode='r')
    stdDev = np.load('%s/specUnc.npy' % directory, mmap_mode='r')
else:
    spectra = np.load('%s/specCube.npy' % directory)
    stdDev = np.load('%s/specUnc.npy' % directory)

h_map = np.load('%s/param.npy' % directory)
vis = np.load('%s/visual.npy' % directory)

h_map[4] = (1./(np.exp(h_map[4]))/60.)

    
global marker
global count
global mask_val
global toggle, toggle2
marker = 1
count = 0
mask_val = 0.005
toggle, toggle2 = 0, 0

global freqs

### determine frequency values that FFT will evaluate
## use frequencies array if exists
if os.path.isfile('%s/frequencies.npy' % directory):
    freqs = np.load('%s/frequencies.npy' % directory)
else:
    if wavelength in [1600, 1700]:
        time_step = 24
    else:
        time_step = 12
    
    freq_size = (spectra.shape[2]*2)+1
    sample_freq = fftpack.fftfreq(freq_size, d=time_step)
    pidxs = np.where(sample_freq > 0)    
    freqs = sample_freq[pidxs]


# assign equal weights to all parts of the curve
df = np.log10(freqs[1:len(freqs)]) - np.log10(freqs[0:len(freqs)-1])
df2 = np.zeros_like(freqs)
df2[0:len(df)] = df
df2[len(df2)-1] = df2[len(df2)-2]
ds0 = df2

 
# create list of figure titles and colorbar names
titles = ['Power Law Slope Coeff.', 'Power Law Index', 'Tail', 
          'Lorentzian Amplitude', 'Lorentzian Location [min]', 
          'Lorentzian Width', 'F-Statistic', 'Averaged Visual Image']
date_title = '%i/%02i/%02i' % (int(date[0:4]),int(date[4:6]),int(date[6:8]))


# create figure with heatmap and spectra side-by-side subplots
fig1 = plt.figure(figsize=(20,10))

ax1 = plt.gca()
ax1 = plt.subplot2grid((30,31),(4, 1), colspan=14, rowspan=25)
ax1.set_xlim(0, h_map.shape[2]-1)
ax1.set_ylim(0, h_map.shape[1]-1)  
ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[1]), 
              y = 1.01, fontsize=17)

param = h_map[1] 
h_min = np.percentile(param,1) 
h_max = np.percentile(param,99)
im, = ([ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, 
                   vmax=h_max, picker=True)])

# design colorbar for heatmaps
global cbar1
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="3%", pad=0.07)
cbar1 = plt.colorbar(im,cax=cax1)
cbar1.ax.tick_params(labelsize=13, pad=3)   

# make toggle buttons to display each parameter's heatmap
axcoeff = plt.axes([0.01, 0.9, 0.05, 0.063])
axindex = plt.axes([0.07, 0.9, 0.05, 0.063])
axtail = plt.axes([0.13, 0.9, 0.05, 0.063])
axlorentz_amp = plt.axes([0.19, 0.9, 0.05, 0.063])
axlorentz_loc = plt.axes([0.25, 0.9, 0.05, 0.063])
axlorentz_wid = plt.axes([0.31, 0.9, 0.05, 0.063])
axfstat = plt.axes([0.37, 0.9, 0.05, 0.063])
axvisual = plt.axes([0.43, 0.9, 0.05, 0.063])
axhist = plt.axes([0.49, 0.9, 0.05, 0.063])
axmask = plt.axes([0.55, 0.9, 0.05, 0.063])
axslider = plt.axes([0.64, 0.915, 0.15, 0.03])
axsaveFig = plt.axes([0.91, 0.9, 0.05, 0.063])

# set up spectra subplot
ax2 = plt.subplot2grid((30,31),(4, 17), colspan=13, rowspan=24)
ax2.loglog()
ax2.set_xlim(10**-4.5, 10**-1.3)
ax2.set_ylim(10**-5, 10**0)  

fig1.canvas.mpl_connect('button_press_event', onclick)

ax2.set_title('Spectra Fit: Pixel ( , )', fontsize=15)
ax2.set_xlabel('Frequency [Hz]', fontsize=fontSz-3, labelpad=5)
ax2.set_ylabel('Power', fontsize=fontSz-3, labelpad=5)
plt.tight_layout()

# add callbacks to each button - linking corresponding action
callback = Index()

bcoeff = Button(axcoeff, 'Coeff.')
bcoeff.on_clicked(callback.coeff)
bindex = Button(axindex, 'Index')
bindex.on_clicked(callback.index)
btail = Button(axtail, 'Tail')
btail.on_clicked(callback.tail)
blorentz_amp = Button(axlorentz_amp, 'Lorentz. Amp')
blorentz_amp.on_clicked(callback.lorentz_amp)
blorentz_loc = Button(axlorentz_loc, 'Lorentz. Loc')
blorentz_loc.on_clicked(callback.lorentz_loc)
blorentz_wid = Button(axlorentz_wid, 'Lorentz. Wid')
blorentz_wid.on_clicked(callback.lorentz_wid)
bfstat = Button(axfstat, 'F-Stat')
bfstat.on_clicked(callback.fstat)
bvisual = Button(axvisual, 'Visual')
bvisual.on_clicked(callback.visual)
bhist = Button(axhist, 'Hist.')
bhist.on_clicked(callback.hist)
bmask = Button(axmask, 'Mask')
bmask.on_clicked(callback.mask)
bsaveFig = Button(axsaveFig, 'Save')
bsaveFig.on_clicked(callback.saveFig)

slid_mask = Slider(axslider, 'Masking', 0.001, 0.1, valinit=0.005)
slid_mask.on_changed(update)

plt.show()