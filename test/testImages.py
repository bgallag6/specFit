import numpy as np
import os
import datetime
from PIL import Image

N = 256 # Number of images
dataDir = './tmp'

if not os.path.exists(dataDir): os.makedirs(dataDir)

###########################################################################
# Create images
imarray = np.zeros((3,3), dtype=np.uint8)
T = float(N)
for i in range(0,N):
    arr = imarray
    t = float(i)
    tmp = np.sin(4.*2*np.pi*t/T)+np.sin(5.*2*np.pi*t/T)+np.sin(6.*2*np.pi*t/T)
    arr[0,0] = np.uint8( 255.0*(tmp/3.+1.)/2. )
    im = Image.fromarray(arr)
    im.save(os.path.join(dataDir,'testImages-%03d.tiff' % i))

###########################################################################
# Read images
cube = np.zeros((3,3,N), dtype=np.uint8)

timestamps = []
exposures = []
for i in range(0,N):
    # Normally we would read timestamp from Exif metadata in file
    timestamp = datetime.datetime.utcfromtimestamp(i).strftime('%Y-%m-%d %H:%M:%S')
    timestamps.append(timestamp)
    
    exposures.append(1.0)

    im = Image.open(os.path.join(dataDir,'testImages-%03d.tiff' % i))
    cube[:,:,i] = np.asarray(im)

cube_avg = np.uint8(np.average(cube,axis=2))

###########################################################################
# Save info for specFit processing
np.save(os.path.join(dataDir,'dataCube.npy'), cube)
np.save(os.path.join(dataDir,'visual.npy'), cube_avg)
np.save(os.path.join(dataDir,'timestamps.npy'), timestamps)
np.save(os.path.join(dataDir,'exposures.npy'),exposures)
