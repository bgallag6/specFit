# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 10:19:24 2018

@author: Brendan
"""

import glob
import numpy as np
import sunpy.cm
from PIL import Image
from sunpy.map import Map
import piexif  # python 2.6+ only?
import PIL.ExifTags
import datetime
from astropy.time import Time

## load FITS image, save as jpg, write timestamp & exposure to EXIF data
def fits2jpg(filename):
    fmap = Map("%s" % filename)
    dmap = fmap.data
    dmap[dmap > 750] = 750
    dmap = dmap//3
    jmap = np.dstack((dmap, dmap, dmap))
    jmap = Image.fromarray(jmap.astype('uint8'))
    timestamp = fmap.date.strftime("%Y:%m:%d %H:%M:%S")
    exposure = (fmap.exposure_time).value
    print(Time(fmap.date).jd, exposure)
    exif_ifd = {piexif.ExifIFD.DateTimeOriginal: u"%s" % timestamp,
            piexif.ExifIFD.ExposureTime: tuple((int(exposure*1e4), int(1e4)))}
    exif_bytes = piexif.dump({"Exif":exif_ifd})
    fnTrim = '%s_%s_%s' % (filename[54:58], filename[45:53], filename[85:-20])
    jmap.save('./jpg2/%s.jpg' % fnTrim, exif=exif_bytes)

## load jpg images, convert to grayscale "matching" original FITS image   
def convertGrayscale(filename):
    im = np.array(Image.open(filename))
    redIm = np.array(im[:,:,0].astype('int16'))
    blueIm = np.array(im[:,:,1].astype('int16'))
    greenIm = np.array(im[:,:,2].astype('int16'))
    grayIm = redIm + blueIm + greenIm
    return grayIm

## extract timestamp & exposure from EXIF data
def readExif(filename):
    im = Image.open(filename)
    exif = {PIL.ExifTags.TAGS[k]: v for k, v in im._getexif().items()
            if k in PIL.ExifTags.TAGS}
    ts = datetime.datetime.strptime(exif["DateTimeOriginal"], fmt)
    timeStamp = Time(ts).jd
    expTime = float(exif["ExposureTime"][0]) / float(exif["ExposureTime"][1])
    return timeStamp, expTime


fmt = "%Y:%m:%d %H:%M:%S"

flist = sorted(glob.glob('./fits2/aia*.fits'))

#for fname in flist:
#    fits2jpg(fname)


