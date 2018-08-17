# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 10:19:24 2018

@author: Brendan
"""

import glob
import sunpy.cm
from PIL import Image

flist = sorted(glob.glob('C:/Users/Brendan/Desktop/specFit/test/Images/20120606/1600/aia*.fits'))

def fits2jpg(filename):
    dmap = Map("%s" % filename).data
    dmap[dmap > 750] = 750
    dmap = dmap//3
    jmap = np.dstack((dmap, dmap, dmap))
    jmap = Image.fromarray(jmap.astype('uint8'))
    fnTrim = '%s_%s_%s' % (filename[54:58], filename[45:53], filename[85:-20])
    jmap.save('C:/Users/Brendan/Desktop/demojpeg/%s.jpg' % fnTrim)
    
def convertGrayscale(filename):
    im = np.array(Image.open(filename))
    redIm = np.array(im[:,:,0].astype('int16'))
    blueIm = np.array(im[:,:,1].astype('int16'))
    greenIm = np.array(im[:,:,2].astype('int16'))
    grayIm = redIm + blueIm + greenIm
    return grayIm

fits2jpg(flist[0])