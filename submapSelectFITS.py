# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 10:47:23 2016

@author: Brendan
"""

import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import glob
from sunpy.map import Map
import astropy.units as u
from astropy.coordinates import SkyCoord
from matplotlib.widgets import Button, TextBox
from sunpy.coordinates import frames


import argparse
parser = argparse.ArgumentParser(description='submapSelect.py')
parser.add_argument('--raw_dir', type=str)

args = parser.parse_args()
raw_dir = args.raw_dir


def onselect(eclick, erelease):
    'eclick and erelease are the press and release events'
    global x1, y1, x2, y2
    global px1, px2, py1, py2
    px1, py1 = eclick.xdata, eclick.ydata
    px2, py2 = erelease.xdata, erelease.ydata
    #print(px1,py1)
  
    sub1 = full_map.pixel_to_world(px1*u.pixel, py1*u.pixel)
    sub2 = full_map.pixel_to_world(px2*u.pixel, py2*u.pixel)
    
    x1 = int(sub1.Tx.value)
    y1 = int(sub1.Ty.value)
    x2 = int(sub2.Tx.value)
    y2 = int(sub2.Ty.value)

    sub_map = full_map.submap(sub1, sub2)
    ax2 = plt.subplot2grid((1,33),(0, 19), colspan=14, rowspan=1, projection=sub_map)
    #ax2 = plt.subplot(122, projection=sub_map)
    sub_map.plot()
    ax2.set_autoscale_on(False)

    t_x1.set_text('Bottom Left: (%ix, %iy)' % (x1, y1))
    t_x2.set_text('Top Right: (%ix, %iy)' % (x2, y2))
    
    plt.draw()
    toggle_selector.RS.update()
  

def saveCoords(event):
        global x1, y1, x2, y2
        print("Saving coordinates: (%i, %i), (%i, %i)." % (x1, y1, x2, y2), flush=True)
        """
        ### placeholder for saving coordinates to config file
        """
        
        
def submit(text):
    global x1, y1, x2, y2
    assert len(list(map(int, text.split()))) == 4, "Please enter 4 coordinates"
    x1, y1, x2, y2 = list(map(int, text.split()))

    sub1 = SkyCoord(x1*u.arcsec, y1*u.arcsec, frame=frames.Helioprojective)
    sub2 = SkyCoord(x2*u.arcsec, y2*u.arcsec, frame=frames.Helioprojective)
    
    sub_map = full_map.submap(sub1, sub2)
    ax2 = plt.subplot2grid((1,33),(0, 19), colspan=14, rowspan=1, projection=sub_map)
    #ax2 = plt.subplot(122, projection=sub_map)
    sub_map.plot()
    ax2.set_autoscale_on(False)
    
    t_x1.set_text('Bottom Left: (%ix, %iy)' % (x1, y1))
    t_x2.set_text('Top Right: (%ix, %iy)' % (x2, y2))

       
def toggle_selector(event):
    pass


def reset(event):
    ax.set_xlim(xlim_default)
    ax.set_ylim(ylim_default)
    plt.draw()
    
"""
############################################
"""

raw_dir = '/Users/bgallagher/Documents/SDO/FITS/20120702/171'
#raw_dir = './images/raw/20120606/1600/fits'

flist = sorted(glob.glob('%s/*.fits' % raw_dir))

nf = len(flist)

mid_file = nf // 2

global full_map
full_map = Map(flist[mid_file])
fm_shape = full_map.data.shape


fig = plt.figure(figsize=(10,6))
plt.subplots_adjust(bottom=0.23)
plt.suptitle('Select a subregion: draw a box or enter coordinates')

ax = plt.subplot2grid((1,33),(0, 0), colspan=14, rowspan=1, projection=full_map)
im = full_map.plot() 
ax.set_autoscale_on(False)

xlim_default = ax.get_xlim()
ylim_default = ax.get_ylim()

toggle_selector.RS = RectangleSelector(ax,onselect,drawtype='box',useblit=True,interactive=True)
plt.connect('key_press_event', toggle_selector)


global x1, y1, x2, y2

x_y1 = full_map.pixel_to_world(0*u.pixel, 0*u.pixel)
x_y2 = full_map.pixel_to_world(fm_shape[1]*u.pixel, fm_shape[0]*u.pixel)

x1 = int(x_y1.Tx.value)
y1 = int(x_y1.Ty.value)
x2 = int(x_y2.Tx.value)
y2 = int(x_y2.Ty.value)

t_x1, = ([fig.text(0.54, 0.085, 'Bottom Left: (%ix, %iy)' % (x1, y1))])
t_x2, = ([fig.text(0.54, 0.055, 'Top Right: (%ix, %iy)' % (x2, y2))])
fig.text(0.242, 0.027, 'Press [Enter] twice after entry', fontsize=8)


axsaveCoords = plt.axes([0.8, 0.04, 0.07, 0.07])
axbox = plt.axes([0.18, 0.05, 0.3, 0.05])
axreset = plt.axes([0.88, 0.04, 0.07, 0.07])


bsaveCoords = Button(axsaveCoords, 'Save')
bsaveCoords.on_clicked(saveCoords)
text_box = TextBox(axbox, 'x1 y1 x2 y2:', initial="")
text_box.on_submit(submit)
breset = Button(axreset, 'Reset')
breset.on_clicked(reset)


plt.draw()
plt.show()