# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 07:34:10 2019

@author: Brendan
"""

import os
import numpy as np

class TestPreProcess():
    
    def test_works(self):
        raw_dir = './images/raw/demo'
        processed_dir = './images/processed/demo'
        Nfiles = 256
        command = ('python preProcessDemo.py --raw_dir {} --processed_dir {} '
                   '--Nfiles {}').format(raw_dir, processed_dir, Nfiles)
        os.system(command)
        assert os.path.exists(os.path.join(processed_dir, 'dataCube.npy'))
    
    def test_cube_dim(self):  
        processed_dir = './images/processed/demo'
        Nfiles = 256
        path = os.path.join(processed_dir, 'dataCube.npy')
        cube_shape = np.load(path).shape
        assert cube_shape[2] == Nfiles

    def test_matching_dim(self):
        processed_dir = './images/processed/demo'
        Nfiles = 256
        ts_size = np.load(os.path.join(processed_dir, 'timestamps.npy')).size
        exp_size = np.load(os.path.join(processed_dir, 'exposures.npy')).size
        assert ts_size == exp_size == Nfiles
    
    def test_periodic(self):
        processed_dir = './images/processed/demo'
        data = np.load(os.path.join(processed_dir, 'dataCube.npy'))[0,0]        
        assert data[0] == data[data.size//2]
        assert data[data.size//8] == data[data.size//8 + data.size//2]
        
        
class TestFftAvg():
    
    def test_works(self):
        processed_dir = './images/processed/demo'
        box_avg_size = 3
        num_segments = 1
        command = ('python fftAvg.py --processed_dir {} --box_avg_size {} '
                   ' --num_segments {}'.format(processed_dir, box_avg_size, num_segments))
        os.system(command)
        assert os.path.exists(os.path.join(processed_dir, 'specCube.npy'))
        
    def test_matching_dim(self):
        processed_dir = './images/processed/demo'
        cube = np.load(os.path.join(processed_dir, 'dataCube.npy'))
        spec = np.load(os.path.join(processed_dir, 'specCube.npy'))
        freqs = np.load(os.path.join(processed_dir, 'frequencies.npy'))
        print("{} (spec) == {} (freqs) == {}/2 (cube)".format(spec.shape[2], freqs.size, (cube.shape[2]//2)))
        assert spec.shape[2] == freqs.size == (cube.shape[2]//2)
        
    def test_freqs(self):
        processed_dir = './images/processed/demo'
        freqs = np.load(os.path.join(processed_dir, 'frequencies.npy'))
        ts = np.load(os.path.join(processed_dir, 'timestamps.npy'))
        tDiff = list(np.diff(ts))
        timeStep = max(tDiff, key=tDiff.count)
        assert (1/freqs[-1]) == (2*timeStep)  # Nyquist frequency
        assert (1/freqs[0]) == (ts.size*timeStep)
    
    def test_power(self):
        processed_dir = './images/processed/demo'
        spec = np.load(os.path.join(processed_dir, 'specCube.npy'))[0,0]
        assert abs(spec.sum() - 1.) < 0.01  # total power very close to 1
       
    