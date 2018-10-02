# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 22:12:20 2018

@author: Brendan
"""

import numpy as np
#import numba  # speed-up ~10%?

# define Power-Law-fitting function (Model M1)
#@numba.jit
def M1(f, a, n, c):
    return a*f**-n + c

# define combined-fitting function (Model M2)
#@numba.jit
def M2(f, a, n, c, p, fp, fw):
    return a*f**-n + c + p*(1./ (1.+((np.log(f)-fp)/fw)**2))

# define Lorentzian-fitting function
#@numba.jit
def m2(f, p, fp, fw):
    return p*(1./ (1.+((np.log(f)-fp)/fw)**2))