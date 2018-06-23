# -*- coding: utf-8 -*-
"""
Created on Thu May 24 08:25:54 2018

@author: Brendan
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import RectangleSelector
import matplotlib.patches as patches
import scipy.signal
from pylab import axvline
import sunpy
from scipy import fftpack
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Button, Slider
from scipy.stats import f as ff
from scipy.stats.stats import pearsonr
import os

    
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
def update(val):
    global mask_val
    mask_val = slid_mask.val
    return mask_val

def plotMap(p):
    param = h_map[p]
    pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
    pNaN = pflat[~np.isnan(pflat)]
    h_min = np.percentile(pNaN,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max = np.percentile(pNaN,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    if p == 4:
        c_map = 'jet_r'
    else:
        c_map = 'jet'
    im = ax1.imshow(param, cmap=c_map, interpolation='nearest', vmin=h_min, vmax=h_max, picker=True)
    ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[p]), y = 1.01, fontsize=17)
    plt.colorbar(im,cax=cax)
    plt.draw()

def plotMask(p):
    global mask_val
    ax1.clear()
    ax1.set_xlim(0, h_map.shape[2]-1)
    ax1.set_ylim(0, h_map.shape[1]-1)  
    
    param = h_map[p]
    pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
    pNaN = pflat[~np.isnan(pflat)]
    h_min = np.percentile(pNaN,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max = np.percentile(pNaN,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    
    # generate p-value heatmap + masked Gaussian component heatmaps
    df1, df2 = 3, 6  # degrees of freedom for model M1, M2
    p_val = ff.sf(h_map[6], df1, df2)
    param_mask = np.copy(param) 
    param_mask[p_val > mask_val] = np.NaN  # mask the Gaussian component arrays with NaNs if above threshold
    
    # determine percentage of region masked 
    count = np.count_nonzero(np.isnan(param_mask))   
    total_pix = p_val.shape[0]*p_val.shape[1]
    mask_percent = ((np.float(count))/total_pix)*100
    
    ax1.set_title(r'%s: %i $\AA$ | %s | $f_{masked}$ = %0.1f%s' % (date_title, wavelength, titles[p], mask_percent, '%'), y = 1.01, fontsize=17)
    if p == 4:
        c_map = 'jet_r'
    else:
        c_map = 'jet'
    im = ax1.imshow(param_mask, cmap=c_map, interpolation='nearest', vmin=h_min, vmax=h_max, picker=True)
    plt.colorbar(im,cax=cax)
    plt.draw()
        
def histMask(p):
    global mask_val
    param = h_map[p]
    
    # generate p-value heatmap + masked Gaussian component heatmaps
    df1, df2 = 3, 6  # degrees of freedom for model M1, M2
    p_val = ff.sf(h_map[6], df1, df2)
    param_mask = np.copy(param) 
    param_mask[p_val > mask_val] = 0.  # mask the Gaussian component arrays with NaNs if above threshold
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
        
    def roll(self, event):  # meh, should probably fix this
        global marker
        marker = 2
        paramA = h_map[0]
        paramn = h_map[1]
        paramC = h_map[2]
        param = (paramC/paramA)**(1./paramn)
        param = np.nan_to_num(param)/60.
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[2]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()      
        
    def lorentz_amp(self, event):
        global marker
        marker = 3
        param = h_map[3]
        pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
        pNaN = pflat[~np.isnan(pflat)]
        h_min = np.percentile(pNaN,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(pNaN,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[3]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
    
    def lorentz_loc(self, event):
        global marker
        marker = 4
        param = h_map[4]
        pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
        pNaN = pflat[~np.isnan(pflat)]
        h_min = np.percentile(pNaN,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(pNaN,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet_r', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[4]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
        
    def lorentz_wid(self, event):
        global marker
        marker = 5
        param = h_map[5]
        pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
        pNaN = pflat[~np.isnan(pflat)]
        h_min = np.percentile(pNaN,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(pNaN,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[5]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
        
    def fstat(self, event):
        global marker
        marker = 6
        param = h_map[6]
        pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
        pNaN = pflat[~np.isnan(pflat)]
        h_min = np.percentile(pNaN,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(pNaN,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='jet', interpolation='none', vmin=h_min, vmax=h_max,  picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[6]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
        
    def visual(self, event):
        param = vis
        h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im = ax1.imshow(param, cmap='sdoaia%i' % wavelength, interpolation='nearest', vmin=h_min, vmax=h_max, picker=True)
        ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[7]), y = 1.01, fontsize=17)
        plt.colorbar(im,cax=cax)
        plt.draw()
    
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
            plt.savefig('%s/%s_%i_%i.pdf' % (outdir, date, wavelength, count), bbox_inches='tight')
        count += 1
        return count
  

# Simple mouse click function to store coordinates
def onclick(event):
    global ix, iy, c
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
        
        # assign equal weights to all parts of the curve
        df = np.log10(f_fit[1:len(f_fit)]) - np.log10(f_fit[0:len(f_fit)-1])
        df2 = np.zeros_like(f_fit)
        df2[0:len(df)] = df
        df2[len(df2)-1] = df2[len(df2)-2]
        ds = df2
        #ds = subcube_StdDev[l][m]
     
        ## fit data to combined power law plus gaussian component model
        
        try:
            # initial guesses for fitting parameters
            M1_low = [-0.002, 0.3, -0.01]
            M1_high = [0.002, 6., 0.01]
            nlfit_l, nlpcov_l = scipy.optimize.curve_fit(PowerLaw, f_fit, s, bounds=(M1_low, M1_high), sigma=ds, method='dogbox')
          
        except RuntimeError:
            print("Error M1 - curve_fit failed")
            pass
        
        except ValueError:
            print("Error M1 - inf/NaN")
            pass
        
        A, n, C = nlfit_l  # unpack fitting parameters
         
       
        try:                                 

            M2_low = [0., 0.3, -0.01, 0.00001, -6.5, 0.05]
            M2_high = [0.002, 6., 0.01, 0.2, -4.6, 0.8]        
                    
            nlfit_gp, nlpcov_gp = scipy.optimize.curve_fit(LorentzPowerBase, f_fit, s, bounds=(M2_low, M2_high), sigma=ds, method='dogbox', max_nfev=3000)
        
        except RuntimeError:
            print("Error M2 - curve_fit failed")
            pass
        
        except ValueError:
            print("Error M2 - inf/NaN")
            pass

        
        A2, n2, C2, P2, fp2, fw2 = nlfit_gp  # unpack fitting parameters
        #print nlfit_gp
               
        try:           
            nlfit_gp2, nlpcov_gp2 = scipy.optimize.curve_fit(LorentzPowerBase, f_fit, s, p0 = [A2, n2, C2, P2, fp2, fw2], bounds=(M2_low, M2_high), sigma=ds, max_nfev=3000)    
           
        except RuntimeError:
            print("Error M2 - curve_fit failed")
            pass
        
        except ValueError:
            print("Error M2 - inf/NaN")
            pass
        
        A22, n22, C22, P22, fp22, fw22 = nlfit_gp2  # unpack fitting parameters     
        #print nlfit_gp2
                       
        # create model functions from fitted parameters    
        m1_fit = PowerLaw(f_fit, A, n, C)    
        lorentz = Lorentz(f_fit,P22,fp22,fw22)
        m2_fit2 = LorentzPowerBase(f_fit, A22,n22,C22,P22,fp22,fw22) 
        m1_fit2 = PowerLaw(f_fit, A22,n22,C2)      
        
        residsM1 = (s - m1_fit)
        chisqrM1 =  ((residsM1/ds)**2).sum() 
        
        residsM22 = (s - m2_fit2)
        chisqrM22 = ((residsM22/ds)**2).sum()
        redchisqrM22 = ((residsM22/ds)**2).sum()/float(f_fit.size-6) 
              
        #f_test = ((chisqrM1-chisqrM2)/(6-3))/((chisqrM2)/(f.size-6))
        f_test2 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(f_fit.size-6))
        
        df1, df2 = 3, 6  # degrees of freedom for model M1, M2
        p_val = ff.sf(f_test2, df1, df2)
        
        amp_scale2 = PowerLaw(np.exp(fp22), A22, n22, C22)  # to extract the gaussian-amplitude scaling factor

        rval = pearsonr(m2_fit2, s)[0]  # calculate r-value correlation coefficient 
        
        rollover = (1. / ((C22 / A22)**(-1. / n22))) / 60.
        
        fwhm = (1. / (np.exp(fp22+fw22) - np.exp(fp22-fw22))) / 60.
        
        ax2.loglog(f_fit, s, 'k', linewidth=1.5)
        
        #ax2.loglog(f, m1_fit, 'r', linewidth=1.3, label='M1')
        ax2.loglog(f_fit, m2_fit2, color='purple', linewidth=1.5, label='M2')
        ax2.loglog(f_fit, m1_fit2, color='green', linewidth=1.5, label='M2: Power Law')
        ax2.loglog(f_fit, lorentz, color='green', linestyle='dashed', linewidth=1.5, label='M2: Lorentz')
        ax2.set_xlabel('Frequency [Hz]', fontsize=font_size-3, labelpad=5)
        ax2.set_ylabel('Power', fontsize=font_size-3, labelpad=5)
        axvline(x=0.00333,color='k',ls='dashed', label='5 minutes')
        axvline(x=0.00555,color='k',ls='dotted', label='3 minutes')
        ax2.set_title('Spectra Fit: Pixel (%ix , %iy)' % (ix, iy), fontsize=font_size-5)
        ax2.text(0.011, 10**-0.5, r'$n$ = {0:0.2f}'.format(n22), fontsize=font_size-2)
        ax2.text(0.011, 10**-0.75, r'$\alpha$ = {0:0.2e}'.format(P22), fontsize=font_size-2)
        ax2.text(0.011, 10**-1.00, r'$\beta$ = {0:0.1f} [min]'.format((1./np.exp(fp22))/60.), fontsize=font_size-2)
        ax2.text(0.0061, 10**-1.25, r'$FWHM$ = {0:0.1f} [min]'.format(fwhm), fontsize=font_size-2)
        #plt.vlines((0.0093),10**-8,10**1, linestyles='dotted', label='3 minutes')
        legend = ax2.legend(loc='lower left', prop={'size':15}, labelspacing=0.35)
        ax2.set_xlim(10**-4.5, 10**-1.3)
        ax2.set_ylim(10**-5, 10**0)   
        ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
        for label in legend.get_lines():
                label.set_linewidth(2.0)  # the legend line width   

    return ix, iy
    
# define combined-fitting function (Model M2)
def LorentzPowerBase(f2, A2, n2, C2, P2, fp2, fw2):
    return A2*f2**-n2 + C2 + P2*(1./ ((np.pi*fw2)*(1.+((np.log(f2)-fp2)/fw2)**2)))

# define Power-Law-fitting function (Model M1)
def PowerLaw(f, A, n, C):
    return A*f**-n + C
        
# define Gaussian-fitting function
def Lorentz(f, P, fp, fw):
    return P*(1./ ((np.pi*fw)*(1.+((np.log(f)-fp)/fw)**2))) 
    

"""
##############################################################################
##############################################################################
"""

#directory = 'S:'
#directory = 'C:/Users/Brendan/Desktop/specFit/test'
directory = './test/validation'
date = '20120606'
wavelength = 1600

global spectra

cube_shape = np.load('%s/Processed/tmp/%s/%i/specCube_shape.npy' % (directory, date, wavelength))
spectra = np.memmap('%s/Processed/tmp/%s/%i/specCube_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))
if wavelength not in [1600, 1700]:
    stddev = np.memmap('%s/Processed/tmp/%s/%i/3x3_stddev_mmap.npy' % (directory, date, wavelength), dtype='float64', mode='r', shape=(cube_shape[0], cube_shape[1], cube_shape[2]))

global marker
global toggle
global toggle2
global count
global mask_val
marker = 1
toggle = 0
toggle2 = 0
count = 0
mask_val = 0.005


### determine frequency values that FFT will evaluate
if wavelength == 1600 or wavelength == 1700:
    time_step = 24
else:
    time_step = 12
freq_size = (cube_shape[2]*2)+1
sample_freq = fftpack.fftfreq(freq_size, d=time_step)
pidxs = np.where(sample_freq > 0)    

plt.rcParams["font.family"] = "Times New Roman"
font_size = 20

if 1:
    
    global f_fit
    
    freqs = sample_freq[pidxs]
    #print(len(freqs))
    f_fit = np.linspace(freqs[0],freqs[len(freqs)-1],int(spectra.shape[2]))   
    
    
    h_map = np.load('%s/Processed/Output/%s/%i/param.npy' % (directory, date, wavelength))
 
    vis = np.load('%s/Processed/Output/%s/%i/visual.npy' % (directory, date, wavelength))

    h_map[4] = (1./(np.exp(h_map[4]))/60.)
    
    # create list of titles and colorbar names for display on the figures
    titles = ['Power Law Slope Coeff.', 'Power Law Index', 'Rollover [min]', 'Lorentzian Amplitude', 'Lorentzian Location [min]', 'Lorentzian Width', 'F-Statistic', 'Averaged Visual Image']
    date_title = '%i/%02i/%02i' % (int(date[0:4]),int(date[4:6]),int(date[6:8]))
    
    # create figure with heatmap and spectra side-by-side subplots
    fig1 = plt.figure(figsize=(20,10))
    ax1 = plt.gca()
    ax1 = plt.subplot2grid((30,31),(4, 1), colspan=14, rowspan=25)

    ax1.set_xlim(0, h_map.shape[2]-1)
    ax1.set_ylim(0, h_map.shape[1]-1)  
    ax1.set_title(r'%s: %i $\AA$ | %s' % (date_title, wavelength, titles[1]), y = 1.01, fontsize=17)
    
    # was getting error "'AxesImage' object is not iterable"
    # - found: "Each element in img needs to be a sequence of artists, not a single artist."
    param = h_map[1]  # set initial heatmap to power law index     
    h_min = np.percentile(param,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max = np.percentile(param,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    im, = ([ax1.imshow(param, cmap='jet', interpolation='nearest', vmin=h_min, vmax=h_max,  picker=True)])
    
    # design colorbar for heatmaps
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = plt.colorbar(im,cax=cax)
    #cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
    cbar.ax.tick_params(labelsize=13, pad=3)   
    
    
    # make toggle buttons to display each parameter's heatmap
    axcoeff = plt.axes([0.01, 0.9, 0.05, 0.063])
    axindex = plt.axes([0.07, 0.9, 0.05, 0.063])
    axroll = plt.axes([0.13, 0.9, 0.05, 0.063])
    axlorentz_amp = plt.axes([0.19, 0.9, 0.05, 0.063])
    axlorentz_loc = plt.axes([0.25, 0.9, 0.05, 0.063])
    axlorentz_wid = plt.axes([0.31, 0.9, 0.05, 0.063])
    axfstat = plt.axes([0.37, 0.9, 0.05, 0.063])
    axvisual = plt.axes([0.43, 0.9, 0.05, 0.063])
    #axscatter = plt.axes([0.49, 0.9, 0.05, 0.063])
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
    ax2.set_xlabel('Frequency [Hz]', fontsize=font_size-3, labelpad=5)
    ax2.set_ylabel('Power', fontsize=font_size-3, labelpad=5)
    plt.tight_layout()
    
    
    # add callbacks to each button - linking corresponding action
    callback = Index()
    
    bcoeff = Button(axcoeff, 'Coeff.')
    bcoeff.on_clicked(callback.coeff)
    bindex = Button(axindex, 'Index')
    bindex.on_clicked(callback.index)
    broll = Button(axroll, 'Rollover')
    broll.on_clicked(callback.roll)
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
    
plt.draw()