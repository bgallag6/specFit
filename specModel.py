# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 22:12:20 2018

@author: Brendan
"""

import numpy as np

# define Power-Law-fitting function (Model M1)
def M1(f, a, n, c):
    return a*f**-n + c

# define combined-fitting function (Model M2)
def M2(f, a, n, c, p, fp, fw):
    return a*f**-n + c + p*(1./ (1.+((np.log(f)-fp)/fw)**2))

# define Lorentzian-fitting function
def Lorentz(f, p, fp, fw):
    return p*(1./ (1.+((np.log(f)-fp)/fw)**2))