# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 22:12:20 2018

@author: Brendan
"""

import numpy as np
#import numba  # speed-up ~10%?

# define power-law-fitting function (Model M1)
#@numba.jit
def M1(f, p0, p1, p2):
    return p0*f**-p1 + p2

# define combined-fitting function (Model M2)
#@numba.jit
def M2(f, p0, p1, p2, p3, p4, p5):
    return p0*f**-p1 + p2 + p3*(1./ (1.+((np.log(f)-p4)/p5)**2))

# define Lorentzian-fitting function
#@numba.jit
def m2(f, p3, p4, p5):
    return p3*(1./ (1.+((np.log(f)-p4)/p5)**2))