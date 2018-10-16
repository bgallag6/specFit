# -*- coding: utf-8 -*-
"""
Created on Thu May 24 08:25:54 2018

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit as Fit
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

    
def ax2setup():
    global title, p_index, p_amp, p_loc, p_wid
    global curveSpec, curveM2, curveM1, curveLorentz
    title, = ([ax2.set_title('Spectra Fit: Pixel ( , )', fontsize=fontSize)])
    ax2.set_xlabel('Frequency [Hz]', fontsize=fontSize, labelpad=5)
    ax2.set_ylabel('Power', fontsize=fontSize, labelpad=5)
    ax2.set_ylim(ylow, yhigh)
    ax2.set_xlim(xlow, xhigh) 
    
    curveSpec, = ax2.loglog(freqs, emptyLine, 'k', linewidth=1.5)
    #curveM1A, = ax2.loglog(freqs, emptyLine, 'r', linewidth=1.3, label='M1')
    curveM2, = ax2.loglog(freqs, emptyLine, c='purple', lw=1.5, label='M2')
    curveM1, = ax2.loglog(freqs, emptyLine, c='g', lw=1.5, label='M2: Power Law')
    curveLorentz, = ax2.loglog(freqs, emptyLine, c='g', ls='--', lw=1.5, label='M2: Lorentz')
    ax2.axvline(x=(1./300.), color='k', ls='dashed', label='5 minutes')
    ax2.axvline(x=(1./180.), color='k', ls='dotted', label='3 minutes')
    
    legend = ax2.legend(loc='lower left', prop={'size':15}, labelspacing=0.35)   
    for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
            
    p_index, = ([ax2.text(0.010, 10**-0.5, '', fontsize=fontSize)])
    p_amp, = ([ax2.text(0.010, 10**-0.75, '', fontsize=fontSize)])
    p_loc, = ([ax2.text(0.010, 10**-1.00, '', fontsize=fontSize)])
    p_wid, = ([ax2.text(0.00635, 10**-1.25, '', fontsize=fontSize)])
    
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
        plotMap(marker)
        return marker       

    def index(self, event):
        global marker
        marker = 1
        plotMap(marker)
        return marker     
        
    def tail(self, event):
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


def specFit(s, ds):
    ## fit data to combined power law plus gaussian component model       
    try:
        m1_param = Fit(M1, freqs, s, p0=M1_guess, bounds=(M1_low, M1_high), 
                      sigma=ds, method='dogbox')[0]
      
    except RuntimeError: print("Error M1 - curve_fit failed")
    except ValueError: print("Error M1 - inf/NaN")
    
    A, n, C = m1_param  # unpack model parameters
     
   
    try:                                                          
        m2_param0 = Fit(M2, freqs, s, p0=M2_guess, bounds=(M2_low, M2_high), 
                       sigma=ds, method='dogbox', max_nfev=3000)[0]
    
    except RuntimeError: print("Error M2 - curve_fit failed")
    except ValueError: print("Error M2 - inf/NaN")

    
    #A2, n2, C2, P2, fp2, fw2 = m2_param0
    #print nlfit_gp
           
    try:   
        m2_param = Fit(M2, freqs, s, p0=m2_param0, bounds=(M2_low, M2_high),
                        sigma=ds, max_nfev=3000)[0]    
       
    except RuntimeError: print("Error M2 - curve_fit failed")
    except ValueError: print("Error M2 - inf/NaN")
    
    A22, n22, C22, P22, fp22, fw22 = m2_param     
    #print m2_param
                   
    # create models from parameters    
    m1_fit = M1(freqs, *m1_param)    
    lorentz = m2(freqs, P22, fp22, fw22)
    m2_fit2 = M2(freqs, *m2_param) 
    m1_fit2 = M1(freqs, A22, n22, C22)      
    
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
    
    #curveM1A.set_ydata(m1_fit)
    curveM2.set_ydata(m2_fit2)
    curveM1.set_ydata(m1_fit2)
    curveLorentz.set_ydata(lorentz)
    
    p_index.set_text(r'$n$ = {0:0.2f}'.format(n22))
    p_amp.set_text(r'$\alpha$ = {0:0.2e}'.format(P22))
    p_loc.set_text(r'$\beta$ = {0:0.1f} [min]'.format(pLoc))
    p_wid.set_text(r'$FWHM$ = {0:0.1f} [min]'.format(fwhm))  


# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy
    ixx, iyy = event.xdata, event.ydata
    if ixx > 1. and iyy > 1.:
        ax2.clear()
        del ax1.collections[:]
        plt.draw()
        ax2setup()
        
        print("pixel: (%ix, %iy)" % (ixx, iyy))
        ix = int(ixx)
        iy = int(iyy)
        
        s = np.array(spectra[iy][ix])
        
        if specVis_fit == True:
            # use 3x3 pixel-box std.dev. or adhoc method for fitting uncertainties
            if spec_unc == 'stddev':
                ds = np.zeros((spectra.shape[2]))
                ds[:] = stdDev[iy][ix][:]
            elif spec_unc == 'adhoc':
                ds = ds0
            specFit(s, ds)
        
        # update subplots
        ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
        
        title.set_text('Spectra Fit: Pixel (%ix , %iy)' % (ix, iy))
        
        curveSpec.set_ydata(s)
                
    return ix, iy

    

"""
##############################################################################
##############################################################################
"""

print("Starting specVis...", flush=True)

plt.rcParams["font.family"] = "Times New Roman"
fontSize = 15


#directory = 'C:/Users/Brendan/Desktop/specFit/images/processed/20120606/1600'
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
specVis_fit = False

global spectra

if mmap_spectra == True:
    # load memory-mapped array as read-only
    spectra = np.load('%s/specCube.npy' % directory, mmap_mode='r')
    stdDev = np.load('%s/specUnc.npy' % directory, mmap_mode='r')
else:
    spectra = np.load('%s/specCube.npy' % directory)
    stdDev = np.load('%s/specUnc.npy' % directory)

vis = np.load('%s/visual.npy' % directory)

haveParam = False
if os.path.isfile('%s/param.npy' % directory):
    h_map = np.load('%s/param.npy' % directory)
    haveParam = True
    
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
date_title = '%i-%02i-%02i' % (int(date[0:4]),int(date[4:6]),int(date[6:8]))


# create figure with heatmap and spectra side-by-side subplots
fig1 = plt.figure(figsize=(18,9))

ax1 = plt.gca()
ax1 = plt.subplot2grid((30,31),(4, 1), colspan=14, rowspan=25)

if haveParam:
    ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[1]), 
                  y = 1.01, fontsize=17)
    h_map[4] = (1./(np.exp(h_map[4]))/60.)
    param = h_map[1] 
    h_min = np.percentile(param,1) 
    h_max = np.percentile(param,99)
    im, = ([ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, 
                       vmax=h_max, picker=True)])
    
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
    slid_mask = Slider(axslider, 'Masking', 0.001, 0.1, valinit=0.005)
    slid_mask.on_changed(update)
    
else:
    ax1.set_title(r'%s: %i $\AA$ | Visual Average' % (date_title, wavelength), 
                  y = 1.01, fontsize=17)
    param = vis
    h_min = np.percentile(param,1) 
    h_max = np.percentile(param,99)
    im, = ([ax1.imshow(param, cmap='sdoaia1600', interpolation='nearest', 
                       vmin=h_min, vmax=h_max, picker=True)])
 
ax1.set_xlim(0, param.shape[1]-1)
ax1.set_ylim(0, param.shape[0]-1)    

# design colorbar for heatmaps
global cbar1
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="3%", pad=0.07)
cbar1 = plt.colorbar(im,cax=cax1)
cbar1.ax.tick_params(labelsize=13, pad=3)   

axsaveFig = plt.axes([0.91, 0.9, 0.05, 0.063])

# add callbacks to each button - linking corresponding action
callback = Index()

bsaveFig = Button(axsaveFig, 'Save')
bsaveFig.on_clicked(callback.saveFig)

fig1.canvas.mpl_connect('button_press_event', onclick)


# set up spectra subplot
ax2 = plt.subplot2grid((30,31),(4, 17), colspan=13, rowspan=24)

xspan = np.log10(freqs[-1]) - np.log10(freqs[0])
xlow = 10**(np.log10(freqs[0]) - (xspan/10))
xhigh = 10**(np.log10(freqs[-1]) + (xspan/10))

yspan = np.log10(np.percentile(spectra, 99.9)) - np.log10(np.percentile(spectra, 0.1))
ylow = 10**(np.log10(np.percentile(spectra, 0.1)) - (yspan/10))
yhigh = 10**(np.log10(np.percentile(spectra, 99.9)) + (yspan/10))

emptyLine = [0 for i in range(len(freqs))]

ax2setup()

plt.tight_layout()

plt.show()