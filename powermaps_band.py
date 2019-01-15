#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 13:38:36 2019

@author: bgallagher
"""

"""
# work on
- why have to click on pixel after changing heatmap?
- When switching back to parameters, don’t necessarily delete heatmap buttons?, also don’t necessarily switch back to visual image
- Maybe include a reset for image (sometimes gets stuck) (take from submapSelectFITS)
- copy any upgrades from htPlot to here:
    - could key-binding be used for anything
    - variable vmin/vmax?
    - change colormap option?
    - self.fig.canvas.draw_idle()
    - have where if press mask and ax2=hist, updates hist (overplot or plot new)
    - maybe don't change hist limits when changing from mask--nomask
    - when have powermaps displayed, when click on pixel - spectra comes up
    - possibly have histogram use subregion?
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit as Fit
import sunpy
import sunpy.cm
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Button, Slider, RadioButtons
from scipy.stats import f as ff
import os
import yaml
from specModel import M1, M2, m2
import glob
from sunpy.map import Map

import argparse
parser = argparse.ArgumentParser(description='specVis.py')
parser.add_argument('--processed_dir', type=str, default='images/processed/demo')
#parser.add_argument('--processed_dir', type=str, default='images/processed/20120606/1600')
parser.add_argument('--raw_dir', type=str, default='images/processed/demo')
args = parser.parse_args()

raw_dir = args.raw_dir
processed_dir = args.processed_dir


def vminVmax(arr):
    # compute 1st and 99th percentile heatmap values (to exclude outliers) 
    flat_arr = arr.flatten()
    no_nan_arr = flat_arr[~np.isnan(flat_arr)]
    min0 = np.percentile(no_nan_arr, 1.)
    max0 = np.percentile(no_nan_arr, 99.)
    return min0, max0


## create plot object that all actions will modify
class cPlot(object):
    
    mask_val = 0.005
    spec_hist = 'spec'
    mask_bool = False
    plot_vers = 0
    marker = 1  # which parameter heatmap is displayed
    
    def __init__(self):
       
        self.fig1 = plt.figure(figsize=(18,10))
        
        # (3) main program options -- buttons in top right
        self.axpre = plt.axes([0.82, 0.93, 0.045, 0.045])
        self.axpost = plt.axes([0.88, 0.93, 0.045, 0.045])
        self.axpower = plt.axes([0.94, 0.93, 0.045, 0.045])
        self.bpre = Button(self.axpre, 'Show\nParameters')
        self.bpre.on_clicked(self.pre)
        self.bpost = Button(self.axpost, 'Show\nTimeseries')
        self.bpost.on_clicked(self.timeSeries)
        self.bpower = Button(self.axpower, 'Show\nPowermaps')
        self.bpower.on_clicked(self.powerMap)
        
        self.fig1.canvas.mpl_connect('button_press_event', onclick)
        self.fig1.canvas.mpl_connect('pick_event', functions.onpick) 
        
        # ax1: displays visual image / heatmaps & location of selected pixel        
        self.ax1 = plt.gca()
        self.ax1 = plt.subplot2grid((31,31),(4, 1), colspan=14, rowspan=16)
        
        h_min, h_max = vminVmax(vis)  # trim bottom/top 1% values (outliers)
        
        if fits:
            c_map = 'sdoaia%i' % wavelength
            title_text = r'%s: %i $\AA$ | Visual Average' % (date, wavelength)
        else:
            c_map = 'viridis'
            title_text = r'Visual Average'
        
            
        self.im = self.ax1.imshow(vis, cmap=c_map, interpolation='nearest', vmin=h_min, 
                       vmax=h_max, picker=True)
        
        self.ax1_title = self.ax1.set_title(title_text, y = 1.01, fontsize=17)
        self.ax1.set_xlim(-0.5, vis.shape[1]-0.5)
        self.ax1.set_ylim(-0.5, vis.shape[0]-0.5) 
        
        # design colorbar for heatmaps
        self.divider = make_axes_locatable(self.ax1)
        self.cax1 = self.divider.append_axes("right", size="3%", pad=0.07)
        self.cbar1 = plt.colorbar(self.im,cax=self.cax1)
        self.cbar1.ax.tick_params(labelsize=13, pad=3)   
        
        self.ax2setup()
    
    def colorBar(self):
        # update heatmap colorbar
        cb1_ticks = np.linspace(self.im.get_clim()[0], self.im.get_clim()[1], 10)
        self.cbar1.set_clim(cb1_ticks[0], cb1_ticks[-1])
        self.cbar1.set_ticks(cb1_ticks) 
        self.cbar1.draw_all()
        self.fig1.canvas.draw_idle()
        
    def colorBar2(self):
        # update power map colorbar
        cb2_ticks = np.linspace(self.im2.get_clim()[0], self.im2.get_clim()[1], 10)
        self.cbar2.set_clim(cb2_ticks[0], cb2_ticks[-1])
        self.cbar2.set_ticks(cb2_ticks) 
        self.cbar2.draw_all()
               
    def ax2setup(self):
        # set up spectra subplot
        global lined, legline, lines
        
        # ax2: displays pixel spectra / parameter histogram
        self.ax2 = plt.subplot2grid((31,31),(4, 17), colspan=13, rowspan=16)
        
        self.ax2_title, = ([self.ax2.set_title('Pixel ( , )', fontsize=fontSize)])
        self.ax2.set_xlabel('Frequency [Hz]', fontsize=fontSize, labelpad=5)
        self.ax2.set_ylabel('Power', fontsize=fontSize, labelpad=5)
        self.ax2.set_ylim(ylow, yhigh)
        self.ax2.set_xlim(xlow, xhigh) 
        
        #import pdb; pdb.set_trace()            
    
        self.curveSpec, = self.ax2.loglog(freqs, emptyCurve, 'k', linewidth=1.5)
        #curveM1A, = ax2.loglog(freqs, emptyCurve, 'r', linewidth=1.3, label='M1')
        self.curveM2, = self.ax2.loglog(freqs, emptyCurve, c='r', lw=1.5, label='M2')
        self.curveM1, = self.ax2.loglog(freqs, emptyCurve, c='g', lw=1.5, label='M2: Power Law')
        self.curveLorentz, = self.ax2.loglog(freqs, emptyCurve, c='g', ls='--', lw=1.5, label='M2: Lorentzian')
        
        self.leg = self.ax2.legend(loc='lower left')
        self.leg.get_frame().set_alpha(0.4)
        
        lines = [self.curveM2, self.curveM1, self.curveLorentz]
        lined = dict()
        
        for legline, origline in zip(self.leg.get_lines(), lines):
            legline.set_picker(5) 
            lined[legline] = origline                
        
    def pre(self, event):
        global l
    
        if cp.plot_vers == 2:
            self.ax3.remove()
            if haveParam:
                for button in self.axbutton:
                    button.remove()
                self.axslider.remove()
                self.axauto.remove()
                self.mask_text.remove()
                self.slid_mask_text.remove()
            #fx.visual()
            self.ax2setup()
        
        if cp.plot_vers == 3:    
            self.ax4.remove()
            self.cbar2.remove()
            self.axslider2.remove()
            if haveParam:
                for button in self.axbutton:
                    button.remove()
                self.axslider.remove()
                self.axauto.remove()
                self.mask_text.remove()
                self.slid_mask_text.remove()
            #fx.visual()
            self.ax2setup()
    
        if cp.plot_vers != 1: 
            cp.plot_vers = 1
            
            if 'ix' not in globals():
                global ix, iy
                ix = spectra.shape[1]//2
                iy = spectra.shape[0]//2
                
            spec = np.array(spectra[iy][ix])
            
            self.ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
            
            self.ax2_title.set_text('Pixel (%ix , %iy): Spectra' % (ix, iy))
            self.curveSpec.set_ydata(spec)
        
            param = fx.paramSliders()
            
            self.curveM2.set_ydata(M2(freqs, *param))
            
            plt.text(0.05, 11.5, "*Using parameters found in: '%s'" % param_dir)
            
            self.axreload = plt.axes([0.83, 0.18, 0.05, 0.05])
            self.axsaveFig = plt.axes([0.83, 0.11, 0.05, 0.05])
            self.axreset = plt.axes([0.83, 0.04, 0.05, 0.05])
            
            self.breload = Button(self.axreload, 'Reload')
            self.breload.on_clicked(fx.reload)
            self.bsaveFig = Button(self.axsaveFig, 'Save')
            self.bsaveFig.on_clicked(fx.saveFig)
            self.breset = Button(self.axreset, 'Reset')
            self.breset.on_clicked(fx.reset)
            
    def timeSeries(self, event):
    
        if cp.plot_vers == 1:
            for slider in fx.axsliders:
                slider.remove()
            
        if cp.plot_vers == 2:  # maybe just: if cp.plot_vers != 2:
            pass
        
        else:  
            if cp.plot_vers == 3:
                self.ax4.remove()
                self.cbar2.remove()
                self.axslider2.remove()
                if haveParam:
                    for button in self.axbutton:
                        button.remove()
                    self.axslider.remove()
                    self.axauto.remove()
                    self.mask_text.remove()
                    self.slid_mask_text.remove()
                #fx.visual()
            
            cp.plot_vers = 2
            
            if 'ix' not in globals():
                global ix, iy
                ix = spectra.shape[1]//2
                iy = spectra.shape[0]//2
            
            self.ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
            
            # ax3: displays pixel timeseries
            self.ax3 = plt.subplot2grid((31,31),(22, 1), colspan=29, rowspan=8)
            self.ax3_title = self.ax3.set_title('Pixel (%ix , %iy): Timeseries' % (ix, iy), fontsize=fontSize)
            self.ax3.set_xlabel('Time', fontsize=fontSize, labelpad=5)
            self.ax3.set_ylabel('Intensity', fontsize=fontSize, labelpad=5)
            self.ax3.set_xlim(timestamps[0]-0.01*t_range, timestamps[-1]+0.01*t_range) 
            
            self.ts, = self.ax3.plot(timestamps, emptyTimeseries, 'k')
            
            self.ax2setup()
       
            timeseries = np.array(imCube[iy+1][ix+1] / exposures)
              
            self.ts.set_ydata(timeseries)  
            self.ax3.set_ylim(timeseries.min()*0.9, timeseries.max()*1.1) 
                
            spec = np.array(spectra[iy][ix])
                
            self.ax2_title.set_text('Pixel (%ix , %iy): Spectra' % (ix, iy))
            self.curveSpec.set_ydata(spec)
            
            if haveParam:
                self.axbutton = []
                
                # make toggle buttons to display each parameter's heatmap
                #axspace = (0.6-(0.01*len(buttons)))/len(buttons)
                for i in range(len(buttons)):
                    self.axbutton.append(plt.axes([0.01+(0.06*i), 0.93, 0.05, 0.045]))
                self.axslider = plt.axes([0.64, 0.935, 0.15, 0.03])
                
                #bFunctions = ['coeff', 'index', 'tail', 'lorentz_amp', 'lorentz_loc', 'lorentz_wid', 'fstat', 'visual', 'hist', 'mask']
                
                self.fnbutton = []
                
                bcount = 0  # use enumerate
                for button in buttons:
                    self.fnbutton.append(Button(self.axbutton[bcount], button))
                    #fnbutton[bcount].on_clicked(eval('fx.%s' % bFunctions[bcount]))
                    self.fnbutton[bcount].on_clicked(fx.ax_loc)
                    bcount += 1
                    
                self.axauto = plt.axes([0.57, 0.93, 0.03, 0.045])
                self.rauto = RadioButtons(self.axauto, ('Off', 'On'))
                self.rauto.on_clicked(fx.maskOnOff)
                self.mask_text = self.fig1.text(0.575, 0.98, 'Mask')    
                
                #slid_mask = Slider(self.axslider, 'P-Value', 0.001, 0.1, valinit=0.005)
                self.slid_mask = Slider(self.axslider, '', 0.001, 0.1, valinit=0.005)
                self.slid_mask.on_changed(functions.update)
                self.slid_mask_text = self.fig1.text(0.70, 0.97, 'P-Value')
                
                self.fig1.canvas.mpl_connect('button_press_event', onclick)   
    
    
    def powerMap(self, event):
        
        if cp.plot_vers != 3:
            if cp.plot_vers != 2:
                self.timeSeries(event)
             
            self.ax3.remove()
            self.ax2.remove()
        
            cp.plot_vers = 3    
            
            # ax4: displays spectra power maps
            self.ax4 = plt.subplot2grid((31,31),(4, 17), colspan=13, rowspan=16)
            self.ax4.set_xlim(0, spectra.shape[1]-1)
            self.ax4.set_ylim(0, spectra.shape[0]-1)  
            self.ax4_title, = ([self.ax4.set_title(r'Period: %0.2f [min]' % 4., y = 1.01, fontsize=17)])
            
            idx_low = (np.abs(freqs - (1./(4*60) - 0.0005))).argmin()
            idx_high = (np.abs(freqs - (1./(4*60) + 0.0005))).argmin()
            
            param = np.copy(spectra[:,:,idx_low:idx_high])  # set initial heatmap to power law index  
            param = np.sum(param, axis=2)
            h_min, h_max = vminVmax(param)

            #self.im2 = self.ax4.imshow(param, cmap='Greys', interpolation='nearest', vmin=h_min, vmax=h_max) 
            self.im2 = self.ax4.imshow(param, cmap='gray', interpolation='nearest', vmin=h_min, vmax=h_max)
            
            # design colorbar for heatmaps
            self.divider2 = make_axes_locatable(self.ax4)
            self.cax2 = self.divider2.append_axes("right", size="3%", pad=0.07)
            self.cbar2 = plt.colorbar(self.im2,cax=self.cax2)
            self.cbar2.ax.tick_params(labelsize=13, pad=3)   
            
            self.axslider2 = plt.axes([0.4, 0.2, 0.3, 0.04])
            #slid_mask = Slider(axslider, 'Frequency', f_fit[0], f_fit[-1], valinit=(1./240))
            self.slid_freqs = Slider(self.axslider2, 'Period', (1./freqs[-1])/60., 50., valinit=4., valfmt='%0.2f')
            self.slid_freqs.on_changed(functions.update3)
            
            self.axslider3 = plt.axes([0.4, 0.13, 0.3, 0.04])
            self.slid_width = Slider(self.axslider3, 'Band Width [mHz]', 0.1, 10., valinit=1., valfmt='%0.2f')
            self.slid_width.on_changed(functions.update4)
        

    

class functions():  
    ind = 0
    
    # ---- params ---- #
    def reload(self, event):
        #also needs to reset -- why does having 'event' here work?
        global M1_low, M1_high, M2_low, M2_high
        with open('specFit_config_test.yaml', 'r') as stream:
            cfg = yaml.load(stream)
        
        M1_low = cfg['M1_low']
        M1_high = cfg['M1_high']
        M2_low = cfg['M2_low']
        M2_high = cfg['M2_high']
        
        for slider in fx.axsliders:
            slider.remove()
        
        param = fx.paramSliders()
        
    def saveFig(self, event):
        print('save params')
        
    def reset(self, event):
        for slider in fx.fnsliders:
            slider.reset()
            
        params = [slider.val for slider in fx.fnsliders] 
    
        cp.curveM2.set_ydata(M2(freqs, *params))
        
    ## make definintion -- used 3 times (be careful as middle use has self.update2 vs functions.update2)
    def paramSliders(self):
        if haveParam:
            param = h_map[iy,ix,:6]
        else:
            param = (np.array(M2_low) + np.array(M2_high)) / 2
        
        self.axsliders = []
        self.fnsliders = []
    
        # make parameter sliders
        for i, M2_label in enumerate(M2_labels):
            self.axsliders.append(plt.axes([0.15, 0.23-(0.04*i), 0.6, 0.02]))
            self.fnsliders.append(Slider(self.axsliders[i], M2_label, M2_low[i], M2_high[i], param[i]))
            self.fnsliders[i].on_changed(functions.update2)
            self.fnsliders[i].valtext.set_text(text[i] % param[i])
        return param
        
    def update2(self):
        params = np.zeros((len(fx.axsliders)))
    
        for i in range(len(fx.axsliders)):
            params[i] = fx.fnsliders[i].val
            fx.fnsliders[i].valtext.set_text(text[i] % params[i])
         
        cp.curveM2.set_ydata(M2(freqs, *params))
        cp.curveM1.set_ydata(M1(freqs, *params[:3]))
        cp.curveLorentz.set_ydata(m2(freqs, *params[3:6]))  
      
    
    # ---- timeseries / power-maps ---- #
    def ax_loc(self, event):
        # execute action based on button clicked
        for i in range(len(cp.axbutton)):
            if event.inaxes == cp.axbutton[i]:
                if i == (len(cp.axbutton)-1):
                    self.hist()
                elif i == (len(cp.axbutton)-2):
                    self.visual()  # maybe assign visual a marker?
                else:
                    cp.marker = i
                    self.plotMap(cp.marker)
        
    def plotMap(self, p):
        param = h_map[:,:,p]
        h_min, h_max = vminVmax(param)

        if p == 4:
            c_map = 'jet_r'
        else:
            c_map = 'jet'
            
        if cp.mask_bool:
            # generate p-value heatmap + masked Lorentzian component heatmaps
            dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
            p_val = ff.sf(h_map[:,:,6], dof1, dof2)
            p_val[np.isnan(p_val)] = 1
            param_mask = np.copy(param) 
            param_mask[p_val > cp.mask_val] = np.NaN
            
            # determine percentage of region masked 
            count = np.count_nonzero(np.isnan(param_mask))   
            total_pix = p_val.shape[0]*p_val.shape[1]
            mask_percent = ((np.float(count))/total_pix)*100
            
            cp.ax1_title.set_text(r'%s | $f_{masked}$ = %0.1f%s' % (titles[p], mask_percent, '%'))
            cp.im.set_data(param_mask)
        
        elif not cp.mask_bool:
            cp.ax1_title.set_text(r'%s' % titles[p])
            cp.im.set_data(param)
            
        cp.im.set_clim(h_min, h_max)
        cp.im.set_cmap(c_map)
        
        cp.colorBar()
    
    """       
    def plotMask(self, p):
        param = h_map[:,:,p]
        h_min, h_max = vminVmax(param)
        
        # generate p-value heatmap + masked Lorentzian component heatmaps
        dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
        p_val = ff.sf(h_map[:,:,6], dof1, dof2)
        p_val[np.isnan(p_val)] = 1
        param_mask = np.copy(param) 
        param_mask[p_val > cp.mask_val] = np.NaN
        
        # determine percentage of region masked 
        count = np.count_nonzero(np.isnan(param_mask))   
        total_pix = p_val.shape[0]*p_val.shape[1]
        mask_percent = ((np.float(count))/total_pix)*100
        
        if p == 4:
            c_map = 'jet_r'
        else:
            c_map = 'jet'
            
        cp.ax1_title.set_text(r'%s | $f_{masked}$ = %0.1f%s' % (titles[p], mask_percent, '%'))
        cp.im.set_data(param_mask)
        cp.im.set_clim(h_min, h_max)
        cp.im.set_cmap(c_map)
        
        cp.colorBar()
    """
    
        
    def visual(self):
        h_min, h_max = vminVmax(vis)

        cp.ax1_title.set_text(r'%s' % titles[7])
        cp.im.set_data(vis)
        cp.im.set_clim(h_min, h_max)
        cp.im.set_cmap('sdoaia%i' % wavelength)
        cp.colorBar()
    
            
    def maskOnOff(self, label):
        label_dict = {'On': True, 'Off': False}
        cp.mask_bool = label_dict[label]
        self.plotMap(cp.marker)
        #if cp.mask_bool:
        #    self.plotMask(cp.marker)
        #elif not cp.mask_bool:
        #    self.plotMap(cp.marker)
     
        
    def update(val):
        cp.mask_val = cp.slid_mask.val
        
    def hist(self):
        if cp.spec_hist != 'hist':
                cp.ax2.set_yscale('linear')
                cp.ax2.set_xscale('linear')
                cp.curveSpec.remove()
                cp.curveM2.remove()
                cp.curveM1.remove()
                cp.curveLorentz.remove()
                cp.spec_hist = 'hist'
                cp.ax2.set_ylabel('Count')
                
        param = h_map[:,:,cp.marker]      
        
        if not cp.mask_bool:
            cp.ax2_title.set_text('Histogram: %s' % titles[cp.marker])
            pflat = param.flatten()
            h_color = 'black'
            
        elif cp.mask_bool:    
            cp.ax2_title.set_text('Histogram: %s | Masked' % titles[cp.marker])
            h_color = 'red'
            
            # ---- generate p-value heatmap + masked Lorentzian component heatmaps
            dof1, dof2 = 3, 6  # degrees of freedom for model M1, M2
            p_val = ff.sf(h_map[:,:,6], dof1, dof2)
            p_val[np.isnan(p_val)] = 1
            param_mask = np.copy(param) 
            param_mask[p_val > cp.mask_val] = 0.
            pflat = param_mask.flatten()
            pflat = pflat[pflat != 0]
                 
        pNaN = pflat[~np.isnan(pflat)]
        y, x, _ = cp.ax2.hist(pNaN, bins=25, edgecolor='black', alpha=0.75, color=h_color)  # need a set_data
        cp.ax2.set_xlabel('%s' % titles[cp.marker])
        cp.ax2.set_xlim(np.percentile(pNaN, 1), np.percentile(pNaN, 99))
        cp.ax2.set_ylim(0, y.max()*1.1)
        cp.leg.set_visible(False)
        plt.draw()    
    
    
    # ---- power-maps ---- #
    def update3(self):
        freq_val = cp.slid_freqs.val
        width_val = cp.slid_width.val
        
        cp.ax4_title.set_text('Period: %0.2f [min]' % freq_val)
        
        freq_val = 1./(freq_val*60)      
        
        freq_val_low = freq_val - (width_val/2000.)  
        freq_val_high = freq_val + (width_val/2000.)
        print('Frequencies: %0.2f - %0.2f [mHz]' % (freq_val_low*1000., freq_val_high*1000.))
        if (freq_val_low < freqs.min()) or (freq_val_high > freqs.max()):
            print('Band outside range of frequencies.')
        if freq_val_low > freq_val_high:
            print('Frequency band inverted.')
        
        idx_low = (np.abs(freqs - freq_val_low)).argmin()
        idx_high = (np.abs(freqs - freq_val_high)).argmin()
        print('Frequency indices: %i-%i' % (idx_low, idx_high))
        
        param = np.copy(spectra[:,:,idx_low:idx_high])
        param = np.sum(param, axis=2)
        h_min, h_max = vminVmax(param)
        
        cp.im2.set_data(param)
        cp.im2.set_clim(h_min, h_max)
        cp.colorBar2()
        
    def update4(self):
        width_val = cp.slid_width.val
        freq_val = 1./(cp.slid_freqs.val*60)
        
        freq_val_low = freq_val - (width_val/2000.)  
        freq_val_high = freq_val + (width_val/2000.)
        print('Frequencies: %0.2f - %0.2f [mHz]' % (freq_val_low*1000., freq_val_high*1000.))
        if (freq_val_low < freqs.min()) or (freq_val_high > freqs.max()):
            print('Band outside range of frequencies.')
        if freq_val_low > freq_val_high:
            print('Frequency band inverted.')
        
        idx_low = (np.abs(freqs - freq_val_low)).argmin()
        idx_high = (np.abs(freqs - freq_val_high)).argmin()
        print('Frequency indices: %i-%i' % (idx_low, idx_high))
        
        param = np.copy(spectra[:,:,idx_low:idx_high])
        param = np.sum(param, axis=2)
        h_min, h_max = vminVmax(param)
        
        cp.im2.set_data(param)
        cp.im2.set_clim(h_min, h_max)
        cp.colorBar2()
    
    
    # ---- toggle visibility of selected lines in legend ---- #
    def onpick(event):
        
        if event.artist in lined:
            legline = event.artist
            origline = lined[legline]
            visb = not origline.get_visible()
            origline.set_visible(visb)
            if visb:
                legline.set_alpha(1.0)
            else:
                legline.set_alpha(0.2)
            cp.fig1.canvas.draw()
        
            

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
    #global fnsliders
    #global axsliders
    ixx, iyy = event.xdata, event.ydata
    
    if event.inaxes == cp.ax1:
        global ix, iy
        del cp.ax1.collections[:]
        plt.draw()
        
        if cp.spec_hist == 'hist':
            cp.ax2setup()
            cp.spec_hist = 'spec'

        #print("location: (%fx, %fy)" % (ixx, iyy))        
        #print("location: (%ix, %iy)" % (ixx, iyy))
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
        cp.ax1.scatter(ix, iy, s=200, marker='x', c='white', linewidth=2.5)
        
        cp.ax2_title.set_text('Pixel (%ix , %iy): Spectra' % (ix, iy))

        cp.curveSpec.set_ydata(s)
    
        timeseries = np.array(imCube[iy+1][ix+1] / exposures)
        
        # update spectra and timeseries
        if cp.plot_vers == 2:
            cp.ts.set_ydata(timeseries)  
            cp.ax3.set_ylim(timeseries.min()*0.9, timeseries.max()*1.1)  
            cp.ax3_title.set_text('Pixel (%ix , %iy): Timeseries' % (ix, iy))
            if haveParam:
                cp.curveM2.set_ydata(M2(freqs, *h_map[iy,ix,:6]))
                cp.curveM1.set_ydata(M1(freqs, *h_map[iy,ix,:3]))
                cp.curveLorentz.set_ydata(m2(freqs, *h_map[iy,ix,3:6]))
        
        # update spectra and model sliders
        if cp.plot_vers == 1:
            
            for slider in fx.axsliders:
                slider.remove()
            
            param = fx.paramSliders()
            
            plt.text(0.05, 11.5, "*Using parameters found in: '%s'" % param_dir)
            
            s = M2(freqs, *param)
            cp.curveM2.set_ydata(s)
            cp.curveM1.set_ydata(M1(freqs, *param[:3]))
            cp.curveLorentz.set_ydata(m2(freqs, *param[3:6]))
            


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
#raw_dir = 'C:/Users/Brendan/Desktop/specFit/images/raw/20120606/1600/fits/'
#processed_dir = '/Users/bgallagher/Documents/SDO/DATA/20120702/1700'
#raw_dir = '/Users/bgallagher/Documents/SDO/FITS/20120702/1700'
processed_dir = '/Users/bgallagher/Documents/External_Transfer/20120308_1700'
raw_dir = '/Users/bgallagher/Documents/External_Transfer/FITS'

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

#buttons = M2_labels + ['F-Stat', 'Visual', 'Hist.', 'Mask']
buttons = M2_labels + ['F-Stat', 'Visual', 'Hist.']

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
          

global freqs
## use frequencies array if exists
if os.path.isfile('%s/frequencies.npy' % processed_dir):
    freqs = np.load('%s/frequencies.npy' % processed_dir)

# assign equal weights to all parts of the curve
df = np.log10(freqs[1:len(freqs)]) - np.log10(freqs[0:len(freqs)-1])
ds0 = np.append(df, df[-1])

def getProperties(filename):
    fmap = Map(filename)
    mapDate = fmap.date.strftime('%Y-%m-%d')
    mapWave = int(fmap.wavelength.value)
    mapXScale = fmap.scale[0].value
    mapYScale = fmap.scale[1].value
    
    return mapDate, mapWave, mapXScale, mapYScale

# create a list of all the fits files
flist = sorted(glob.glob(raw_dir+'/*'))

if flist[0].find('.fits') != -1:
    date, wavelength, xscale, yscale = getProperties(flist[0])
    fits = True

xspan = np.log10(freqs[-1]) - np.log10(freqs[0])
xlow = 10**(np.log10(freqs[0]) - (xspan/10))
xhigh = 10**(np.log10(freqs[-1]) + (xspan/10))


yspan = np.log10(np.percentile(spectra, 99.9)) - np.log10(np.percentile(spectra, 0.1))
ylow = 10**(np.log10(np.percentile(spectra, 0.1)) - (yspan/10))
yhigh = 10**(np.log10(np.percentile(spectra, 99.9)) + (yspan/10))

emptyCurve = [0 for i in range(len(freqs))]
emptyTimeseries = [-1 for i in range(len(timestamps))]


#ax1.set_xticklabels(labels2int(ax1.get_xticks()))
#ax1.set_yticklabels(labels2int(ax1.get_yticks()))
    

#import pdb; pdb.set_trace()            
"""
if param.size <= 100:
    ax1.set_xticks(-0.5+np.arange(0, param.shape[0], 1), minor=True)
    ax1.set_xticks(np.arange(0, param.shape[0], 1))
    ax1.set_yticks(-0.5+np.arange(0, param.shape[0], 1), minor=True)
    ax1.set_yticks(np.arange(0, param.shape[0], 1))

    ax1.grid(which='minor')
"""    
#ax1.set_xlim(-0.5, param.shape[1]-0.5)
#ax1.set_ylim(-0.5, param.shape[0]-0.5)     

t_range = timestamps[-1]-timestamps[0]

fx = functions()
cp = cPlot()

plt.tight_layout()

plt.show()