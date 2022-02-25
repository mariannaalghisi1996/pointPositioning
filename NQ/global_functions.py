# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 11:43:41 2022

@author: maria
"""

import numpy as np

def thirdInterpolationFunction(z1, z2, z3, z4, x):
    if abs(2*x) < 10**(-10):
        zx = z2
    else:
        delta = 2*x - 1
        g1 = z3 + z2
        g2 = z3 - z2
        g3 = z4 + z1
        g4 = (z4 - z1)/3
        
        a0 = 9*g1 - g3
        a1 = 9*g2 - g4
        a2 = g3 - g1
        a3 = g4 - g2
        
        zx = (1/16)*(a0 + a1*delta + a2*(delta**2) + a3*(delta**3))
    return zx

def epst(x, y, z, w):
    ret = x*np.exp((w-y)/z)/(1 + np.exp((w-y)/z))**2
    return ret
