from PIL import Image
import numpy as np
import os
import sys

T = 256 # Number of images
tmpdir = './tmp'

if not os.path.exists(tmpdir): os.makedirs(tmpdir)

imarray = np.zeros((200,200), dtype=np.uint8)
cube = np.zeros((200,200,T), dtype=np.uint8)
for i in range(0,T):
    arr = imarray+np.uint8(i)
    arr[0,0] = np.sin(2*np.pi*i/(4.*T))+np.sin(2*np.pi*i/(5.*T))+np.sin(2*np.pi*i/(6*T))
    cube[:,:,2] = arr
    im = Image.fromarray(arr)
    im.save(os.path.join(tmpdir,'testImages-%03d.tiff' % i))

sys.path.insert(0,'../')
from fftAvg import fftAvg


