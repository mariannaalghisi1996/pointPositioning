# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 15:04:14 2022

@author: maria
"""

from NQ.constants import *
from NQ.global_functions import *
from NQ.modelParams import *

import math
import numpy as np

'''
Electron density computation
'''
def bottomElectronDensity(h, p):
    A1, A2, A3, hmF2, hmF1, hmE, B2bot, B1top, B1bot, BEtop, BEbot = p.get('A1'), p.get('A2'), p.get('A3'), p.get('hmF2'), p.get('hmF1'), p.get('hmE'), p.get('B2bot'), p.get('B1top'), p.get('B1bot'), p.get('BEtop'), p.get('BEbot')
    # Relevant parameters selection:
    if h > hmE:
        BE = BEtop
    else:
        BE = BEbot
    
    if  h > hmF1:
        BF1 = B1top
    else:
        BF1 = B1bot
    
    # Exponential arguments for each layer
    if h < 100:
        halpha = 100
    else:
        halpha = h
    
    alpha1 = (halpha - hmF2)/B2bot
    alpha2 = ((halpha - hmF1)/BF1)*np.exp(10/(1 + abs(halpha - hmF2)))
    alpha3 = ((halpha - hmE)/BE)*np.exp(10/(1 + abs(halpha - hmF2)))
    
    if abs(alpha1) <= 25:
        s1 = A1*np.exp(alpha1)/(1 + np.exp(alpha1))**2
    else:
        s1 = 0
    
    if abs(alpha2) <= 25:
        s2 = A2*np.exp(alpha2)/(1 + np.exp(alpha2))**2
    else:
        s2 = 0
    
    if abs(alpha3) <= 25:
        s3 = A3*np.exp(alpha3)/(1 + np.exp(alpha3))**2
    else:
        s3 = 0
    
    if h >= 100:
        N = (s1 + s2 + s3)*10**11
    else:
        # Compute corrective parameters
        
        if abs(alpha1) <= 25:
            ds1 = (1/B2bot)*(1 - np.exp(alpha1))/(1 + np.exp(alpha1))
        else:
            ds1 = 0
        if abs(alpha2) <= 25:
            ds2 = (1/BF1)*(1 - np.exp(alpha2))/(1 + np.exp(alpha2))
        else:
            ds2 = 0
        if abs(alpha3) <= 25:
            ds3 = (1/BE)*(1 - np.exp(alpha3))/(1 + np.exp(alpha3))
        else:
            ds3 = 0
        
        # Compute Chapman parameters
        
        BC = 1 - 10*(s1*ds1 + s2*ds2 + s3*ds3)/(s1 + s2 + s3)
        z = (h - 100)/10
        
        N = (s1 + s2 + s3)*np.exp(1 - BC*z - np.exp(-z))*10**11
    
    return N

def topElectronDensity(h, p):
    NmF2, hmF2, H0 = p.get('NmF2'), p.get('hmF2'), p.get('H0')
    # Constant parameters
    g = 0.125
    r = 100
    
    delta_h = h - hmF2
    z = delta_h/(H0*(1 + (r*g*delta_h)/(r*H0 + g*delta_h)))
    #print(z)
    ea = np.exp(z)
    #print(ea)
    
    if ea > 10**11:
        N = (4*NmF2/ea)*10**11
    else:
        N = (4*NmF2*ea/(1 + ea)**2)*10**11
    
    return N

def getN(h, params):
    hmF2 = params.get('hmF2')
    
    if h<= hmF2:
        N = bottomElectronDensity(h, params)
    else:
        N = topElectronDensity(h, params)
    return N

