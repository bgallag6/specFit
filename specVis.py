# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 21:05:07 2018

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
import glob
from sunpy.map import Map

import argparse
parser = argparse.ArgumentParser(description='specVis.py')
#parser.add_argument('--processed_dir', type=str, default='images/processed/demo')
parser.add_argument('--processed_dir', type=str, default='images/processed/20120606/1600')
args = parser.parse_args()

processed_dir = args.processed_dir

def ax2setup():
    global title, curveSpec, ax2
    global lined, legline, lines
    # set up spectra subplot
    ax2 = plt.subplot2grid((30,31),(4, 17), colspan=13, rowspan=16)
    
    title, = ([ax2.set_title('Pixel ( , )', fontsize=fontSize)])
    ax2.set_xlabel('Frequency [Hz]', fontsize=fontSize, labelpad=5)
    ax2.set_ylabel('Power', fontsize=fontSize, labelpad=5)
    ax2.set_ylim(ylow, yhigh)
    ax2.set_xlim(xlow, xhigh) 
    
    curveSpec, = ax2.loglog(freqs, emptyLine, 'k', linewidth=1.5)
    #curveM1A, = ax2.loglog(freqs, emptyLine, 'r', linewidth=1.3, label='M1')
    if haveParam:
        global curveM2, curveM1, curveLorentz
        curveM2, = ax2.loglog(freqs, emptyLine, c='r', lw=1.5, label='M2')
        curveM1, = ax2.loglog(freqs, emptyLine, c='g', lw=1.5, label='M2: Power Law')
        curveLorentz, = ax2.loglog(freqs, emptyLine, c='g', ls='--', lw=1.5, label='M2: Lorentzian')
        
        leg = ax2.legend(loc='lower left')
        leg.get_frame().set_alpha(0.4)
        
        lines = [curveM2, curveM1, curveLorentz]
        lined = dict()
        
        for legline, origline in zip(leg.get_lines(), lines):
            legline.set_picker(5) 
            lined[legline] = origline
        

def onpick(event):
    # on the pick event, find the orig line corresponding to the
    # legend proxy line, and toggle the visibility
    #print(event.artist)
    if event.artist in lined:
        legline = event.artist
        origline = lined[legline]
        visb = not origline.get_visible()
        origline.set_visible(visb)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if visb:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig1.canvas.draw()

 

def update2(val):
    params = np.zeros((len(axsliders)))

    for i in range(len(axsliders)):
        params[i] = fnsliders[i].val
        fnsliders[i].valtext.set_text(text[i] % params[i])
        """
        if i == 0:
            params[i] = 10**fnsliders[i].val
        else:
            params[i] = fnsliders[i].val
        if i == 4:
            fnsliders[i].valtext.set_text(text[i] % (1./(np.exp(params[4])*60.)))
        else:
            fnsliders[i].valtext.set_text(text[i] % params[i])
        """
       
    s = M2(freqs, *params)
    
    #l.set_ydata(s)
    curveM2.set_ydata(s)
    curveM1.set_ydata(M1(freqs, *params[:3]))
    curveLorentz.set_ydata(m2(freqs, *params[3:6]))

def reset(event):
    for slider in fnsliders:
        slider.reset()
    c_m2.set_ydata(emptyLine)
    c_l2.set_ydata(emptyLine)
    params = [slider.val for slider in fnsliders]
    s = M2(freqs, *params)
    
    #l.set_ydata(s)
    curveM2.set_ydata(s)
    
def update(val):
    global mask_val
    mask_val = slid_mask.val
    return mask_val

def update3(val):
    global mask_val
    global im2
    cbar2.remove()
    im2.remove()
    #mask_val = np.log(slid_mask.val)
    #slid_mask.valtext.set_text(mask_val)
    mask_val = slid_mask2.val
    mask_val = 1./(mask_val*60)
    ax4.clear()
    ax4.set_xlim(0, h_map.shape[1]-1)
    ax4.set_ylim(0, h_map.shape[0]-1) 
    title4.set_text('Period: %0.2f' % mask_val)
    
    idx = (np.abs(freqs - mask_val)).argmin()
    
    #param = np.zeros((spectra.shape[0], spectra.shape[1]))
    param = np.copy(spectra[:,:,idx])
    
    pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
    pNaN = pflat[~np.isnan(pflat)]
    h_min = np.percentile(pNaN,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
    h_max = np.percentile(pNaN,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
    
    #ax1.set_title(r'%s: %i $\AA$ | %s | $f_{masked}$ = %0.1f%s' % (date_title, wavelength, titles[p], mask_percent, '%'), y = 1.01, fontsize=17)
    im2 = ax4.imshow(param, cmap='Greys', interpolation='nearest', vmin=h_min, vmax=h_max)
    colorBar2()
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

def colorBar2():
    global cax2
    global cbar2
    # design colorbar for heatmaps
    divider2 = make_axes_locatable(ax4)
    cax2 = divider2.append_axes("right", size="3%", pad=0.07)
    cbar2 = plt.colorbar(im2,cax=cax2)
    #cbar.set_label('%s' % cbar_labels[1], size=15, labelpad=10)
    cbar2.ax.tick_params(labelsize=13, pad=3)   
    plt.colorbar(im2,cax=cax2)
    plt.draw()
    

def plotMap(p):
    global cbar1
    global im
    cbar1.remove()
    im.remove()
    param = h_map[:,:,p]
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
    ax1.set_title(r'%s' % titles[p], y = 1.01, fontsize=17)
    colorBar()
    
def plotMask(p):
    global mask_val
    global cbar1
    global im
    cbar1.remove()
    im.remove()  
    
    param = h_map[:,:,p]
    pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
    pNaN = pflat[~np.isnan(pflat)]
    h_min = np.percentile(pNaN,1)
    h_max = np.percentile(pNaN,99)
    
    # generate p-value heatmap + masked Lorentzian component heatmaps
    dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
    p_val = ff.sf(h_map[:,:,6], dof1, dof2)
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
    ax1.set_title(r'%s | $f_{masked}$ = %0.1f%s' % (titles[p], mask_percent, '%'), 
                  y=1.01, fontsize=17)
    colorBar()
        
def histMask(p):
    global mask_val
    global hist0

    param = h_map[:,:,p]
    
    # generate p-value heatmap + masked Lorentzian component heatmaps
    dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
    p_val = ff.sf(h_map[:,:,6], dof1, dof2)
    param_mask = np.copy(param) 
    param_mask[p_val > mask_val] = 0.
    param1d = np.reshape(param_mask, (param_mask.shape[0]*param_mask.shape[1]))
    pmask = param1d[param1d != 0]
    pmask = pmask[~np.isnan(pmask)]
    title.set_text('Histogram: %s | Masked' % titles[marker])
    hist0 = ax2.hist(pmask, bins=25, edgecolor='black')
    
def visual():
    global cbar1
    global im
    cbar1.remove()
    im.remove()
    param = vis
    h_min = np.percentile(param,1)
    h_max = np.percentile(param,99)
    im = ax1.imshow(param, cmap='sdoaia1600', interpolation='nearest', 
                    vmin=h_min, vmax=h_max, picker=True)
    ax1.set_title(r'%s' % titles[7], y = 1.01, fontsize=17)
    colorBar()

def hist():
    global marker
    global toggle2
    global hist0
    #ax2.clear()
    if toggle2 == 0:
        title.set_text('Histogram: %s' % titles[marker])
        pflat = np.reshape(h_map[:,:,marker], (h_map[:,:,marker].shape[0]*h_map[:,:,marker].shape[1]))
        pNaN = pflat[~np.isnan(pflat)]
        hist0 = ax2.hist(pNaN, bins=25, edgecolor='black')
        ax2.set_xlim(0.3,4)
        ax2.set_ylim(0,1000)
    elif toggle2 == 1:
        histMask(marker)
    plt.draw()
  
def mask():
    global toggle2
    global marker
    if toggle2 == 0:
        toggle2 = 1
        plotMask(marker)
    elif toggle2 == 1:
        toggle2 = 0
        plotMap(marker)              
    return toggle2


def setPre():
    global toggle3
    global l, c_m2, c_l2
    global axsaveFig, axreset, axreload, bsaveFig, breset, breload
    global axsliders, fnsliders

    if toggle3 == 2:
        ax3.remove()
        if haveParam:
            for button in axbutton:
                button.remove()
            axslider.remove()
        visual()
    
    if toggle3 == 3:    
        ax4.remove()
        cbar2.remove()
        axslider2.remove()
        if haveParam:
            for button in axbutton:
                button.remove()
            axslider.remove()
        visual()

    if toggle3 != 1: 
        toggle3 = 1
 
        emptyLine = [0 for i in range(len(freqs))]
        
        if 'ix' not in globals():
            global ix, iy
            ix = spectra.shape[1]//2
            iy = spectra.shape[0]//2
            
        ##have so that if no param file, then maybe load middle of param bounds
        spec = np.array(spectra[iy][ix])
            
        title.set_text('Pixel (%ix , %iy): Spectra' % (ix, iy))
        
        curveSpec.set_ydata(spec)
        
        ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
        
        axsliders = []        
        fnsliders = []
        
        if haveParam:
            param = h_map[iy,ix,:6]
        else:
            global curveM2
            param = (np.array(M2_low) + np.array(M2_high)) / 2
            curveM2, = ax2.loglog(freqs, emptyLine, c='r', lw=1.5, label='M2')
        
        #s = M2(freqs, *h_map[iy,ix,:6])
        s = M2(freqs, *param)
        
        # make parameter sliders
        for i, M2_label in enumerate(M2_labels):
            axsliders.append(plt.axes([0.15, 0.23-(0.04*i), 0.6, 0.02]))
            fnsliders.append(Slider(axsliders[i], M2_label, M2_low[i], M2_high[i], param[i]))
            fnsliders[i].on_changed(update2)
            fnsliders[i].valtext.set_text(text[i] % param[i])
            
        curveM2.set_ydata(s)
        #l, = ax2.loglog(freqs, s, lw=1.5, color='red')
        c_m2, = ax2.loglog(freqs, emptyLine, 'b', linewidth=1.3, label='M2 - Lorentz')
        c_l2, = ax2.loglog(freqs, emptyLine, 'b--', linewidth=1.3, label='Lorentz')
        
        plt.text(0.05, 11.5, "*Using parameters found in: '%s'" % param_dir)
         
        axreload = plt.axes([0.83, 0.18, 0.05, 0.05])
        axsaveFig = plt.axes([0.83, 0.11, 0.05, 0.05])
        axreset = plt.axes([0.83, 0.04, 0.05, 0.05])
        
        # add callbacks to each button - linking corresponding action
        callback = Index()
        
        breload = Button(axreload, 'Reload')
        breload.on_clicked(callback.reload)
        bsaveFig = Button(axsaveFig, 'Save')
        bsaveFig.on_clicked(callback.saveFig)
        breset = Button(axreset, 'Reset')
        breset.on_clicked(reset)

def setPost():
    global toggle3, ts, ax3, title3
    global axbutton, axslider, fnbutton, slid_mask

    if toggle3 == 1:
        for slider in axsliders:
            slider.remove()

        #l.remove()
        axsaveFig.remove()
        axreset.remove()
    if toggle3 == 2:
        pass
    else:  
        if toggle3 == 3:
            ax4.remove()
            cbar2.remove()
            axslider2.remove()
            if haveParam:
                for button in axbutton:
                    button.remove()
                axslider.remove()
            visual()
        
        toggle3 = 2
        
        if 'ix' not in globals():
            global ix, iy
            ix = spectra.shape[1]//2
            iy = spectra.shape[0]//2
        
        ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
            
        ax3 = plt.subplot2grid((30,31),(21, 1), colspan=29, rowspan=8)
        title3, = ([ax3.set_title('Pixel (%ix , %iy): Timeseries' % (ix, iy), fontsize=fontSize)])
        ax3.set_xlabel('Time', fontsize=fontSize, labelpad=5)
        ax3.set_ylabel('Intensity', fontsize=fontSize, labelpad=5)
        ax3.set_xlim(timestamps[0]-0.01*t_range, timestamps[-1]+0.01*t_range) 
        #ax3.set_ylim(0,1)
        
        emptyLine2 = [-1 for i in range(len(timestamps))]
        
        ts, = ax3.plot(timestamps, emptyLine2, 'k')
        
        ax2setup()
   
        timeseries = np.array(imCube[iy+1][ix+1] / exposures)
          
        ts.set_ydata(timeseries)  
        ax3.set_ylim(timeseries.min()*0.9, timeseries.max()*1.1) 
            
        spec = np.array(spectra[iy][ix])
            
        title.set_text('Pixel (%ix , %iy): Spectra' % (ix, iy))
        
        curveSpec.set_ydata(spec)
        
        if haveParam:
            axbutton = []
            
            # make toggle buttons to display each parameter's heatmap
            #axspace = (0.6-(0.01*len(buttons)))/len(buttons)
            for i in range(len(buttons)):
                axbutton.append(plt.axes([0.01+(0.06*i), 0.9, 0.05, 0.063]))
            axslider = plt.axes([0.64, 0.915, 0.15, 0.03])
            
            # add callbacks to each button - linking corresponding action
            callback = Index()
            
            #bFunctions = ['coeff', 'index', 'tail', 'lorentz_amp', 'lorentz_loc', 'lorentz_wid', 'fstat', 'visual', 'hist', 'mask']
            
            fnbutton = []
            
            bcount = 0
            for button in buttons:
                fnbutton.append(Button(axbutton[bcount], button))
                #fnbutton[bcount].on_clicked(eval('callback.%s' % bFunctions[bcount]))
                fnbutton[bcount].on_clicked(callback.ax_loc)
                bcount += 1
                
            slid_mask = Slider(axslider, 'Masking', 0.001, 0.1, valinit=0.005)
            slid_mask.on_changed(update)
            
            fig1.canvas.mpl_connect('button_press_event', onclick)   


def setPower():
    global axslider2, slid_mask2, ax4
    global im2, toggle3, title4
    
    if toggle3 != 3:
        if toggle3 != 2:
            setPost()
         
        ax3.remove()
        ax2.remove()
    
        toggle3 = 3    
        
        ax4 = plt.subplot2grid((30,31),(4, 17), colspan=13, rowspan=16)
        ax4.set_xlim(0, h_map.shape[1]-1)
        ax4.set_ylim(0, h_map.shape[0]-1)  
        title4, = ([ax4.set_title(r'Period: %0.2f [min]' % 4., y = 1.01, fontsize=17)])
        
        idx = (np.abs(freqs - 1./(4*60))).argmin()
        param = np.copy(spectra[:,:,idx])  # set initial heatmap to power law index     
        pflat = np.reshape(param, (param.shape[0]*param.shape[1]))
        pNaN = pflat[~np.isnan(pflat)]
        h_min = np.percentile(pNaN,1)  # set heatmap vmin to 1% of data (could lower to 0.5% or 0.1%)
        h_max = np.percentile(pNaN,99)  # set heatmap vmax to 99% of data (could up to 99.5% or 99.9%)
        im2, = ([ax4.imshow(param, cmap='Greys', interpolation='nearest', vmin=h_min, vmax=h_max)])  
        
        global cbar2
        # design colorbar for heatmaps
        divider2 = make_axes_locatable(ax4)
        cax2 = divider2.append_axes("right", size="3%", pad=0.07)
        cbar2 = plt.colorbar(im2,cax=cax2)
        cbar2.ax.tick_params(labelsize=13, pad=3)   
        
        axslider2 = plt.axes([0.4, 0.2, 0.3, 0.04])
        
        #slid_mask = Slider(axslider, 'Frequency', f_fit[0], f_fit[-1], valinit=(1./240))
        #slid_mask = Slider(axslider, 'Period', (1./f_fit[-1])/60., (1./f_fit[0])/60., valinit=4., valfmt='%0.2f')
        slid_mask2 = Slider(axslider2, 'Period', (1./freqs[-1])/60., 50., valinit=4., valfmt='%0.2f')
        slid_mask2.on_changed(update3)
        
    
class Index(object):
    ind = 0
       
    def ax_loc(self, event):
        for i in range(len(axbutton)):
            if event.inaxes == axbutton[i]:
                if i == (len(axbutton)-1):
                    mask()
                elif i == (len(axbutton)-2):
                    hist()
                elif i == (len(axbutton)-3):
                    visual()
                else:
                    global marker
                    marker = i
                    plotMap(marker)
                    return marker       

    def pre(self, event):
        setPre()
        
    def post(self, event):
        setPost()
        
    def power(self, event):
        setPower()
        
    def saveFig(self, event):
        print('save params')
        
    def reload(self, event):
        global fnsliders, axsliders 
        global M1_low, M1_high, M2_low, M2_high
        with open('specFit_config_test.yaml', 'r') as stream:
            cfg = yaml.load(stream)
        
        M1_low = cfg['M1_low']
        M1_high = cfg['M1_high']
        M2_low = cfg['M2_low']
        M2_high = cfg['M2_high']
        
        for slider in axsliders:
            slider.remove()
            
        axsliders = []
        fnsliders = []
        
        if haveParam:
            param = h_map[iy,ix,:6]
        else:
            param = (np.array(M2_low) + np.array(M2_high)) / 2
        
        # make parameter sliders
        for i, M2_label in enumerate(M2_labels):
            axsliders.append(plt.axes([0.15, 0.23-(0.04*i), 0.6, 0.02]))
            fnsliders.append(Slider(axsliders[i], M2_label, M2_low[i], M2_high[i], param[i]))
            fnsliders[i].on_changed(update2)
            fnsliders[i].valtext.set_text(text[i] % param[i])
        
        

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
    
    """
    residsM1 = (s - m1_fit)
    chisqrM1 =  ((residsM1/ds)**2).sum() 
    
    residsM22 = (s - m2_fit2)
    chisqrM22 = ((residsM22/ds)**2).sum()
    redchisqrM22 = ((residsM22/ds)**2).sum()/float(freqs.size-6) 
          
    f_test2 = ((chisqrM1-chisqrM22)/(6-3))/((chisqrM22)/(freqs.size-6))
    
    dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
    #p_val = ff.sf(f_test2, dof1, dof2)
    
    # extract the lorentzian amplitude scaling factor
    amp_scale2 = M1(np.exp(fp22), A22, n22, C22)  
    """
    
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


# select pixel in ax1
def onclick(event):
    #global ix, iy
    global fnsliders
    global axsliders
    ixx, iyy = event.xdata, event.ydata
    
    if event.inaxes == ax1:
        global ix, iy
        del ax1.collections[:]
        plt.draw()
        
        print("location: (%ix, %iy)" % (ixx, iyy))
        ix = int(round(ixx))
        iy = int(round(iyy))
        
        s = np.array(spectra[iy][ix])
        
        if specVis_fit == True:
            # use 3x3 pixel-box std.dev. or adhoc method for fitting uncertainties
            if spec_unc == 'stddev':
                ds = np.array(stdDev[iy][ix])
            elif spec_unc == 'adhoc':
                ds = ds0
            specFit(s, ds)
        
        # update subplots
        ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
        
        title.set_text('Pixel (%ix , %iy): Spectra' % (ix, iy))
        
        curveSpec.set_ydata(s)
    
        timeseries = np.array(imCube[iy+1][ix+1] / exposures)
          
        if toggle3 == 2:
            ts.set_ydata(timeseries)  
            ax3.set_ylim(timeseries.min()*0.9, timeseries.max()*1.1)  
            title3.set_text('Pixel (%ix , %iy): Timeseries' % (ix, iy))
            if haveParam:
                curveM2.set_ydata(M2(freqs, *h_map[iy,ix,:6]))
                curveM1.set_ydata(M1(freqs, *h_map[iy,ix,:3]))
                curveLorentz.set_ydata(m2(freqs, *h_map[iy,ix,3:6]))
            
        if toggle3 == 1:
            
            for slider in axsliders:
                slider.remove()
                
            axsliders = []
            fnsliders = []
            
            if haveParam:
                param = h_map[iy,ix,:6]
            else:
                param = (np.array(M2_low) + np.array(M2_high)) / 2
            
            # make parameter sliders
            for i, M2_label in enumerate(M2_labels):
                axsliders.append(plt.axes([0.15, 0.23-(0.04*i), 0.6, 0.02]))
                fnsliders.append(Slider(axsliders[i], M2_label, M2_low[i], M2_high[i], param[i]))
                fnsliders[i].on_changed(update2)
                fnsliders[i].valtext.set_text(text[i] % param[i])
                
            plt.text(0.05, 11.5, "*Using parameters found in: '%s'" % param_dir)
            
            s = M2(freqs, *param)
            #l.set_ydata(s)
            curveM2.set_ydata(s)
            curveM1.set_ydata(M1(freqs, *param[:3]))
            curveLorentz.set_ydata(m2(freqs, *param[3:6]))
            


def labels2int(labels):
    '''Remove tick labels with fraction.'''
    newlabels = []
    for i in range(0,len(labels)):
        if int(labels[i]) != labels[i]:
            newlabels.append('')
        else:
            newlabels.append(str(int(labels[i])))
    return newlabels
    

"""
##############################################################################
##############################################################################
"""

#processed_dir = 'C:/Users/Brendan/Desktop/specFit/images/processed/20120606/1600'

#print("Starting specVis...", flush=True)

plt.rcParams["font.family"] = "Times New Roman"
#plt.rcParams["font.size"] = 15
fontSize = 15

with open('specFit_config_test.yaml', 'r') as stream:
    cfg = yaml.load(stream)

mmap_spectra = cfg['mmap_spectra']
M1_low = cfg['M1_low']
M1_high = cfg['M1_high']
M2_low = cfg['M2_low']
M2_high = cfg['M2_high']
spec_unc = cfg['spec_unc']
M1_guess = cfg['M1_guess']
M2_guess = cfg['M2_guess']
M2_labels = cfg['M2_labels']
map_titles = cfg['M2_titles']
#processed_dir = cfg['specVis_dir']
specVis_fit = False
fits = False
text = ['%0.2e', '%0.2f', '%0.2e', '%0.2e', '%0.2f', '%0.2f']

titles = map_titles + ['F-Statistic', 'Averaged Visual Image']

buttons = M2_labels + ['F-Stat', 'Visual', 'Hist.', 'Mask']

global spectra

if mmap_spectra == True:
    # load memory-mapped array as read-only
    spectra = np.load('%s/specCube.npy' % processed_dir, mmap_mode='r')
    stdDev = np.load('%s/specUnc.npy' % processed_dir, mmap_mode='r')
    imCube = np.load('%s/dataCube.npy' % processed_dir, mmap_mode='r')
else:
    spectra = np.load('%s/specCube.npy' % processed_dir)
    stdDev = np.load('%s/specUnc.npy' % processed_dir)
    imCube = np.load('%s/dataCube.npy' % processed_dir)
 
vis = np.load('%s/visual.npy' % processed_dir)

timestamps = np.load('%s/timestamps.npy' % processed_dir)
exposures = np.load('%s/exposures.npy' % processed_dir)

haveParam = False
if os.path.isfile('%s/param.npy' % processed_dir):
    h_map = np.load('%s/param.npy' % processed_dir)
    #h_map[:,:,4] = (1./(np.exp(h_map[:,:,4]))/60.)
    haveParam = True
    param_dir = './%s/param.npy' % processed_dir
else:
    param_dir = './specFit_config.yaml'
           
    
global marker
global count
global mask_val
global toggle, toggle2, toggle3
marker = 1
count = 0
mask_val = 0.005
toggle, toggle2, toggle3 = 0, 0, 0

global freqs

### determine frequency values that FFT will evaluate
## use frequencies array if exists
if os.path.isfile('%s/frequencies.npy' % processed_dir):
    freqs = np.load('%s/frequencies.npy' % processed_dir)

# assign equal weights to all parts of the curve
df = np.log10(freqs[1:len(freqs)]) - np.log10(freqs[0:len(freqs)-1])
df2 = np.zeros_like(freqs)
df2[0:len(df)] = df
df2[len(df2)-1] = df2[len(df2)-2]
ds0 = df2


# create figure with heatmap and spectra side-by-side subplots
fig1 = plt.figure(figsize=(18,9))

ax1 = plt.gca()
ax1 = plt.subplot2grid((30,31),(4, 1), colspan=14, rowspan=16)


def getProperties(filename):
    fmap = Map(filename)
    mapDate = fmap.date.strftime('%Y-%m-%d')
    mapWave = int(fmap.wavelength.value)
    mapXScale = fmap.scale[0].value
    mapYScale = fmap.scale[1].value
    
    return mapDate, mapWave, mapXScale, mapYScale

raw_dir = 'C:/Users/Brendan/Desktop/specFit/images/raw/20120606/1600/fits'

# create a list of all the fits files
flist = sorted(glob.glob('%s/aia*.fits' % raw_dir))

if flist[0].find('.fits') != -1:
    date, wavelength, xscale, yscale = getProperties(flist[0])
    fits = True


if fits:
    ax1.set_title(r'%s: %i $\AA$ | Visual Average' % (date, wavelength), y = 1.01, fontsize=17)
    param = vis 
    #param = h_map[:,:,1]
    h_min = np.percentile(param,1) 
    h_max = np.percentile(param,99)
    im, = ([ax1.imshow(param, cmap='sdoaia%i' % wavelength, interpolation='nearest', vmin=h_min, 
                       vmax=h_max, picker=True)])
    
else:
    ax1.set_title(r'Visual Average', y = 1.01, fontsize=17)
    param = vis
    #h_min = np.percentile(param,1) 
    #h_max = np.percentile(param,99)

    h_min = 0
    h_max = 127

    print(param)
    param[0,0] = 10
    param[0,1] = 20
    param[0,2] = 30
    print(param)

    im, = ([ax1.imshow(param, cmap=plt.get_cmap('viridis', 32), 
                       vmin=h_min, vmax=h_max, picker=True)])

    ax1.set_xticklabels(labels2int(ax1.get_xticks()))
    ax1.set_yticklabels(labels2int(ax1.get_yticks()))
    
    
ax1.set_xlim(0, param.shape[1]-0.5)
ax1.set_ylim(0, param.shape[0]-0.5)     

# design colorbar for heatmaps
global cbar1
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="3%", pad=0.07)
cbar1 = plt.colorbar(im,cax=cax1)
cbar1.ax.tick_params(labelsize=13, pad=3)   

# add callbacks to each button - linking corresponding action
callback = Index()

axpre = plt.axes([0.82, 0.91, 0.045, 0.045])
axpost = plt.axes([0.88, 0.91, 0.045, 0.045])
axpower = plt.axes([0.94, 0.91, 0.045, 0.045])
bpre = Button(axpre, 'Show\nParameters')
bpre.on_clicked(callback.pre)
bpost = Button(axpost, 'Show\nTimeseries')
bpost.on_clicked(callback.post)
bpower = Button(axpower, 'Show\nPowermaps')
bpower.on_clicked(callback.power)

fig1.canvas.mpl_connect('button_press_event', onclick)

fig1.canvas.mpl_connect('pick_event', onpick)   

xspan = np.log10(freqs[-1]) - np.log10(freqs[0])
xlow = 10**(np.log10(freqs[0]) - (xspan/10))
xhigh = 10**(np.log10(freqs[-1]) + (xspan/10))

yspan = np.log10(np.percentile(spectra, 99.9)) - np.log10(np.percentile(spectra, 0.1))
ylow = 10**(np.log10(np.percentile(spectra, 0.1)) - (yspan/10))
yhigh = 10**(np.log10(np.percentile(spectra, 99.9)) + (yspan/10))

emptyLine = [0 for i in range(len(freqs))]

ax2setup()


"""
# set up spectra subplot
ax2 = plt.subplot2grid((30,31),(4, 17), colspan=13, rowspan=16)

title, = ([ax2.set_title('Pixel ( , )', fontsize=fontSize)])
ax2.set_xlabel('Frequency [Hz]', fontsize=fontSize, labelpad=5)
ax2.set_ylabel('Power', fontsize=fontSize, labelpad=5)
ax2.set_ylim(ylow, yhigh)
ax2.set_xlim(xlow, xhigh) 

curveSpec, = ax2.loglog(freqs, emptyLine, 'k', linewidth=1.5)
#curveM1A, = ax2.loglog(freqs, emptyLine, 'r', linewidth=1.3, label='M1')
if haveParam:
    curveM2, = ax2.loglog(freqs, emptyLine, c='r', lw=1.5, label='M2')
if specVis_fit:
    curveM2, = ax2.loglog(freqs, emptyLine, c='purple', lw=1.5, label='M2')
    curveM1, = ax2.loglog(freqs, emptyLine, c='g', lw=1.5, label='M2: Power Law')
    curveLorentz, = ax2.loglog(freqs, emptyLine, c='g', ls='--', lw=1.5, label='M2: Lorentz')
    #ax2.axvline(x=(1./300.), color='k', ls='dashed', label='5 minutes')
    #ax2.axvline(x=(1./180.), color='k', ls='dotted', label='3 minutes')    
    legend = ax2.legend(loc='lower left', prop={'size':15}, labelspacing=0.35)   
    for label in legend.get_lines():
            label.set_linewidth(2.0)  # the legend line width
        
    p_index, = ([ax2.text(0.010, 10**-0.5, '', fontsize=fontSize)])
    p_amp, = ([ax2.text(0.010, 10**-0.75, '', fontsize=fontSize)])
    p_loc, = ([ax2.text(0.010, 10**-1.00, '', fontsize=fontSize)])
    p_wid, = ([ax2.text(0.00635, 10**-1.25, '', fontsize=fontSize)])
"""

t_range = timestamps[-1]-timestamps[0]

plt.tight_layout()

plt.show()