# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 10:47:23 2016

@author: Brendan
"""

from __future__ import print_function
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import glob
from sunpy.map import Map
import astropy.units as u
from astropy.coordinates import SkyCoord
from matplotlib.widgets import Button, TextBox
from sunpy.coordinates import frames


class Index(object):
    
    def saveFig(self, event):
        global x1, y1, x2, y2
        print("Saving coordinates: (%i, %i), (%i, %i)." % (x1, y1, x2, y2), flush=True)
        """
        ### placeholder for saving coordinates to config file
        """

def onselect(eclick, erelease):
    'eclick and erelease are the press and release events'
    global x1, y1, x2, y2
    px1, py1 = eclick.xdata, eclick.ydata
    px2, py2 = erelease.xdata, erelease.ydata
  
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


"""
############################################
"""

#raw_dir = 'S:/FITS/20180211/1700'
raw_dir = './images/raw/20120606/1600/fits'

flist = sorted(glob.glob('%s/*.fits' % raw_dir))

nf = len(flist)

mid_file = nf // 2


full_map = Map(flist[mid_file])
fm_shape = full_map.data.shape


fig = plt.figure(figsize=(10,6))
plt.subplots_adjust(bottom=0.23)
plt.suptitle('Select a subregion: draw a box or enter coordinates')

#ax = plt.subplot(121, projection=full_map)
ax = plt.subplot2grid((1,33),(0, 0), colspan=14, rowspan=1, projection=full_map)
im = full_map.plot() 
ax.set_autoscale_on(False)

toggle_selector.RS=RectangleSelector(ax,onselect,drawtype='box',useblit=True,interactive=True)
plt.connect('key_press_event', toggle_selector)


"""
### maybe have full region plotted initially
ax2 = plt.subplot(122, projection=m1)
m1.plot()
ax2.set_autoscale_on(False)
ax2.set_title('Subregion | Middle File')
"""

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


axsaveFig = plt.axes([0.8, 0.04, 0.07, 0.07])
axbox = plt.axes([0.18, 0.05, 0.3, 0.05])

# add callbacks to each button - linking corresponding action
callback = Index()

bsaveFig = Button(axsaveFig, 'Save')
bsaveFig.on_clicked(callback.saveFig)
text_box = TextBox(axbox, 'x1 y1 x2 y2:', initial="")
text_box.on_submit(submit)


plt.draw()
plt.show()