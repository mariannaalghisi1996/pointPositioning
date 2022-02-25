# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 11:20:07 2022

@author: maria
"""

from NQ.CCIR_MoDIP.ccir_fm3 import Fm3
from NQ.CCIR_MoDIP.ccir_f2 import F2
from NQ.CCIR_MoDIP.modip import stModip
from NQ.constants import DR, RD, Zenith_0
from NQ.global_functions import thirdInterpolationFunction, epst

import math
import numpy as np
import georinex as gr

def UTtime(time):
    UT = time.hour + time.minute/60 + time.second/3600 #UT in hours and decimals
    return UT

def LTtime(time, lon):
    LT = UTtime(time) + lon/15
    return LT

def getMODIP(lat, lon):
    '''
    Parameters
    ----------
    lat : latitude of the receiver in degrees.
    lon : longitude of the receiver in degrees.

    Returns
    -------
    MODIP : Modified Dip Latitude at the location of the user receiver in degrees.

    '''
    # Extreme cases
    if (abs(lat) == 90):
        MODIP = lat
    # For all other latitudes we follow this routine
    else:
        l = round((lon + 180)/10) - 2

        if l < 0:
            l = l + 36
        elif l > 33:
            l = l - 36
        
        a = (lat + 90)/5 +1
        x = a - round(a)
        i = round(a) - 2
        
        Z = []
        
        for k in range(4):
            z_k = []
            for j in range(4):
                z_jk = stModip[i+j][l+k]
                z_k.append(z_jk)
            Z.append(z_k)
        
        z_k = []
        for k in range(4):
            z1 = Z[0][k]
            z2 = Z[1][k]
            z3 = Z[2][k]
            z4 = Z[3][k]
            
            zx_k = thirdInterpolationFunction(z1, z2, z3, z4, x)
            z_k.append(zx_k)
        
        b = (lon + 180)/10
        y = b - round(b)
        
        MODIP = thirdInterpolationFunction(z_k[0], z_k[1], z_k[2], z_k[3], y)
        
        return MODIP

def effIonisationLevelAz(lat, lon, nav_path):
    '''
    Parameters
    ----------
    lat : latitude of the reciver in degrees.
    lon : longitude of the receiver in degrees.
    nav_path : path to navigation rinex file, used to retreve Az params (iono correction params).

    Returns
    -------
    Az : Effective Ionization level Az --> index that represent the daily solar activity 
         according to the receiver location and .

    '''
    MODIP = getMODIP(lat, lon)
    
    hdr = gr.rinexheader(nav_path)
    Gparams = (hdr.get('IONOSPHERIC CORR')).get('GAL')
    
    if Gparams[0] == Gparams[1] == Gparams[2] == 0:
        Az = 63.7
    else:
        Az = Gparams[0] + Gparams[1]*MODIP + Gparams[2]*MODIP**2
    
    if (Az > 0 and Az<400):
        return Az
    else:
        raise TypeError('Invalid Az - Effective ionization level')

def effSunspotNumberAzR(lat, lon, nav_path):
    Az_R = np.sqrt(167273 + (effIonisationLevelAz(lat, lon, nav_path)-63.7)*1123.6) - 408.99
    return Az_R

def solarZenith(time, lat, lon):
    '''
    

    Parameters
    ----------
    time : YY/mm/dd HH:MM:SS.
    lat : [degrees].
    lon : [degrees].

    Returns
    -------
    Z : Solar Zenith angle for a given location in degrees.

    '''
    # Computation of the solar declination:
    mth = time.month
    d_y = 30.5*mth - 15
    # Time in days
    t = d_y + (18-UTtime(time))/24
    # Argument:
    a_m = (0.9856*t - 3.289)*DR
    a_l = a_m + (1.916*np.sin(a_m) + 0.020*np.sin(2*a_m) + 282.634)*DR
    # Sin and cos of solar declination
    sin_Dsun = 0.39782*np.sin(a_l)
    cos_Dsun = np.sqrt(1 - sin_Dsun**2)
    # Solar zenith angle computation:
    cos_Z = np.sin(lat*DR)*sin_Dsun + np.cos(lat*DR)*cos_Dsun*np.cos(math.pi/12*(12-LTtime(time, lon)))
    Z = RD*np.arctan2(np.sqrt(1-cos_Z**2), cos_Z) #Risultato in degree
    
    return Z

def effSolarZenith(time, lat, lon):
    Z = solarZenith(time, lat, lon)
    Z_eff = (Z+(90-0.24*np.exp(20-0.2*Z))*np.exp(12*(Z-Zenith_0)))/(1+np.exp(12*(Z-Zenith_0)))
    
    return Z_eff

def ElayerParams(lat, lon, time, nav_path):
    '''
    Parameters
    ----------
    lat : [degrees].
    lon : [degrees].
    time : YY/mm/dd HH:MM:SS.
    nav_path : path to navigation rinex file.

    Returns
    -------
    foE: E layer critical frequency [MHz]
    NmE: E layer maximum density [10^11/m^3]
    '''
    Az = effIonisationLevelAz(lat, lon, nav_path)
    Z_eff = effSolarZenith(time, lat, lon)
    
    # Definition of the seasonal parameter
    mth = time.month
    winter = [1, 2, 11, 12]
    mid = [3, 4, 8, 10]
    summer = [5, 6, 7, 8]

    if mth in winter:
        seas = -1
    elif mth in mid:
        seas = 0
    elif mth in summer:
        seas = 1
    
    ee = np.exp(0.3*lat)
    seasp = seas*(ee-1)/(ee+1)
    
    foE = np.sqrt(((1.112 - 0.019*seasp)**2)*np.sqrt(Az)*(np.cos(Z_eff*DR))**0.6+0.49)
    NmE = 0.124*foE**2
    
    return [foE, NmE]

def F1layerParams(lat, lon, time, nav_path):
    '''
    Parameters
    ----------
    lat : [degrees].
    lon : [degrees].
    time : YY/mm/dd HH:MM:SS.
    nav_path : path to navigation rinex file.

    Returns
    -------
    foF1 = F1 layer critical frequency [MHz]
    NmF1 = F1 layer maximum density NmF1 [10^11/m^3]

    '''
    foE = ElayerParams(lat, lon, time, nav_path)[0]

    if foE >= 2.0:
        foF1 = 1.4*foE
    else:
        foF1 = 0
    # Correction of the parameter
    if foF1 < 10**(-6):
        foF1 = 0
        
    if (foF1 <= 0 and foE > 2):
        NmF1 = 0.124*(foE + 0.5)**2
    else:
        NmF1 = 0.124*foF1**2

    return [foF1, NmF1]

def F2layerParams(lat, lon, time, nav_path):
    '''
    Parameters
    ----------
    lat : [degrees].
    lon : [degrees].
    time : YY/mm/dd HH:MM:SS.
    nav_path : path to navigation rinex file.

    Returns
    -------
    foF2 = F2 layer critical frequency [MHz]
    NmF2 = F2 layer maximum density NmF1 [10^11/m^3]
    M(3000)

    '''
    mth = time.month
    Az_R = effSunspotNumberAzR(lat, lon, nav_path)
    MODIP = getMODIP(lat, lon)
    
    #Compute AF2, the array of interpolated coefficients for foF2 and Am3, the array of interpolated coefficients for M(3000)F2
    
    F2_m = F2.get(mth)
    fm3_m = Fm3.get(mth)
    F2_1 = F2_m[0]
    F2_2 = F2_m[1]
    
    af2 = []
    for i in range(76):
        new_row = []
        for j in range(13):
            x = float(F2_1[i][j])*(1-Az_R/100)+float(F2_2[i][j])*(Az_R/100)
            new_row.append(x)
        af2.append(new_row)
    
    fm3_1 = fm3_m[0]
    fm3_2 = fm3_m[1]
    
    Am3 = []
    for i in range(49):
        new_row = []
        for j in range(9):
            x = float(fm3_1[i][j])*(1-Az_R/100)+float(fm3_2[i][j])*(Az_R/100)
            new_row.append(x)
        Am3.append(new_row)
    
    # Computation of Fourier time series for foF2 and M(3000)
    
    T = (15*UTtime(time) - 180)*DR
    CF2 = []
    for i in range(76):
        somm = 0
        for j in range(6):
            somm = af2[i][2*j]*np.sin(j*T) + af2[i][2*j+1]*np.cos(j*T)
        somm = somm + af2[i][0]
        CF2.append(somm)
    
    Cm3 = []
    for i in range(49):
        somm = 0
        for j in range(4):
            somm = Am3[i][2*j]*np.sin(j*T) + Am3[i][2*j+1]*np.cos(j*T)
        somm = somm + Am3[i][0]
        Cm3.append(somm)
    
    '''foF2 and M(3000)F2 computation'''
    
    # MODIP coefficients
    m = [1]
    for k in range(2,13):
        m_k = np.sin(MODIP*DR)**(k-1)
        m.append(m_k)
    
    # Latitude and Longitude coefficients
    p, s, c = [], [], []
    for n in range(2, 10):
        p.append(np.cos(lat*DR)**(n-1))
        s.append(np.sin((n-1)*lon*DR))
        c.append(np.cos((n-1)*lon*DR))
    
    #foF2 computation
    
    foF2_1 = 0
    for i in range(12):
        foF2_1 = foF2_1 + CF2[i]*m[i]
    
    Q = [12,12,9,5,2,1,1,1,1]
    K = [-Q[0]]
    for i in range(1,10):
        K.append(K[i-1]+2*Q[i-1])
    
    foF2_c = [foF2_1]    
    for n in range(1,8):
        foF2_n = 0
        for k in range(Q[n]):
            foF2_n = foF2_n + (CF2[K[n]+2*k-1]*c[n] + CF2[K[n]+2*k]*s[n])*m[k]*p[n]
        foF2_c.append(foF2_n)
        
    foF2 = 0
    for i in foF2_c:
        foF2 = foF2 + i
    
    #Compute M(3000)/F2
    
    M_3000_F2_1 = 0
    for k in range(7): 
        M_3000_F2_1 = M_3000_F2_1 + Cm3[k]*m[k]
        
    R = [7,8,6,3,2,1,1]
    
    H = [-R[0]]
    for n in range(1, 7):
        H.append(H[n-1] + 2*R[n-1])
    
    M_3000_F2_c = [M_3000_F2_1]
    for n in range(1,7):
        M_3000_F2_n = 0
        for k in range(R[n]):
            M_3000_F2_n = M_3000_F2_n + (Cm3[H[n]+2*k-1]*c[n]+Cm3[H[n]+2*k]*s[n])*m[k]*p[n]
        M_3000_F2_c.append(M_3000_F2_n)
    
    M_3000 = 0
    
    for i in M_3000_F2_c:
        M_3000 = M_3000 + i
    
    # Compute NmF2
    
    NmF2 = 0.124*foF2**2
    
    return [foF2, NmF2, M_3000]

def getLayerMaxDensityHeight(lat, lon, time, nav_path):
    '''
    Parameters
    ----------
    lat : [degrees].
    lon : [degrees].
    time : YY/mm/dd HH:MM:SS.
    nav_path : path to navigation rinex file.

    Returns
    -------
    The function compute the maximum density height for each layer: hmE, hmF1, hmF2 [km].    

    '''
    # E layer
    foE = ElayerParams(lat, lon, time, nav_path)[0]
    hmE = 120
    
    # F2 layer
    P = F2layerParams(lat, lon, time, nav_path)
    foF2, M = P[0], P[2]
    
    # F1 layer
    ro = (
      (foF2/foE)*np.exp(20*(foF2/foE-1.75)) + 1.75
      )/(np.exp(20*(foF2/foE-1.75))+1)

    if foE < 10**(-30):
        deltaM = -0.012
    else:
        deltaM = 0.253/(ro - 1.215) - 0.012
    
    hmF2 = 1490*M*np.sqrt((0.0196*M**2+1)/(1.2967*M**2-1))/(M+deltaM)-176
    hmF1 = (hmE + hmF2)/2
    
    return [hmE, hmF1, hmF2]

def thicknessParams(lat, lon, time, nav_path):
    '''
    Parameters
    ----------
    lat : [degrees].
    lon : [degrees].
    time : YY/mm/dd HH:MM:SS.
    nav_path : path to navigation rinex file.

    Returns
    -------
    The function returns thickness parameters B2bot, B1top, B1bot, BEtop, BEbot [km]
    
    '''
    PF2 = F2layerParams(lat, lon, time, nav_path)
    foF2, NmF2, M = PF2[0],PF2[1], PF2[2]
    
    PH = getLayerMaxDensityHeight(lat, lon, time, nav_path)
    hmE, hmF1, hmF2 = PH[0], PH[1], PH[2]
    
    B2bot = 0.385*NmF2/(0.01*np.exp(-3.467 + 0.857*np.log(foF2**2)+2.02*np.log(M)))
    B1top = 0.3*(hmF2 - hmF1)
    B1bot = 0.5*(hmF1 - hmE)
    BEtop = max(B1bot, 7)
    BEbot = 5
    
    return [B2bot, B1top, B1bot, BEtop, BEbot]

def getAmplitude(lat, lon, time, nav_path):
    '''
    Parameters
    ----------
    lat : [degrees].
    lon : [degrees].
    time : YY/mm/dd HH:MM:SS.
    nav_path : path to navigation rinex file.

    Returns
    -------
    A1 = F2 layer amplitude [10**11/m**3]
    A2 = F1 layer amplitude [10**11/m**3]
    A3 = E layer amplitude [10**11/m**3]

    '''
    # Required inputs
    E_l = ElayerParams(lat, lon, time, nav_path)
    F1_l = F1layerParams(lat, lon, time, nav_path)
    F2_l = F2layerParams(lat, lon, time, nav_path)
    hm_l = getLayerMaxDensityHeight(lat, lon, time, nav_path)
    B_l = thicknessParams(lat, lon, time, nav_path)
    
    NmF2 = F2_l[1]
    foF1 = F1_l[0]
    NmF1 = F1_l[1]
    NmE = E_l[1]
    hmF2 = hm_l[2]
    hmF1 = hm_l[1]
    hmE = hm_l[0]
    B2bot = B_l[0]
    B1bot = B_l[2]
    BEtop = B_l[3]
    
    A1 = 4*NmF2
    
    if foF1 < 0.5:
        A2 = 0
        A3 = 4*(NmE - epst(A1, hmF2, B2bot, hmE))
    else:
        A3a = 4*NmE
        for iterate in range(5):
            A2a = 4.0*(NmF1 - epst(A1, hmF2, B2bot, hmF1) - epst(A3a, hmE, BEtop, hmF1))
            A2a = (A2a*np.exp(A2a - 0.8*NmF1) + 0.8*NmF1)/(1 + np.exp(A2a - 0.8*NmF1))
            A3a = 4*(NmE - epst(A2a, hmF1, B1bot, hmE) - epst(A1, hmF2, B2bot, hmE))
        A2 = A2a
        A3 = (A3a*np.exp(60*(A3a-0.005))+0.05)/(1 + np.exp(60*(A3a-0.005)))
    
    return [A1, A2, A3]
    
def getH0(lat, lon, time, nav_path):
    mth = time.month
    Az_R = effSunspotNumberAzR(lat, lon, nav_path)
    NmF2 = F2layerParams(lat, lon, time, nav_path)[1]
    hmF2 = getLayerMaxDensityHeight(lat, lon, time, nav_path)[2]
    B2bot = thicknessParams(lat, lon, time, nav_path)[0]
    
    # Shape paramether k
    if mth in [4, 5, 6, 7, 8, 9]:
        ka = 6.705 - 0.014*Az_R - 0.008*hmF2
    else:
        ka = -7.77 + 0.097*(hmF2/B2bot)**2 + 0.153*NmF2
    kb = (ka*np.exp(ka - 2) + 2)/(1 + np.exp(ka - 2))
    k = (8*np.exp(kb - 8) + kb)/(1 + np.exp(kb - 8))
    
    # Compute the topside thickness parameter H0 
    Ha = k*B2bot
    x = (Ha - 150)/100
    v = (0.041163*x - 0.183981)*x + 1.424472
    H0 = Ha/v
    
    return H0
    
def getModelParams(lat, lon, time, nav_path):
    E_l = ElayerParams(lat, lon, time, nav_path)
    F1_l = F1layerParams(lat, lon, time, nav_path)
    F2_l = F2layerParams(lat, lon, time, nav_path)
    hm_l = getLayerMaxDensityHeight(lat, lon, time, nav_path)
    B_l = thicknessParams(lat, lon, time, nav_path)
    A_l = getAmplitude(lat, lon, time, nav_path)
    params = {
        'foE': E_l[0], 
        'NmE': E_l[1],
        'foF1': F1_l[0],
        'NmF1': F1_l[1],
        'foF2': F2_l[0], 
        'NmF2': F2_l[1], 
        'M_3000': F2_l[2],
        'hmE': hm_l[0], 
        'hmF1': hm_l[1], 
        'hmF2': hm_l[2],
        'B2bot': B_l[0], 
        'B1top': B_l[1], 
        'B1bot': B_l[2], 
        'BEtop': B_l[3], 
        'BEbot': B_l[4],
        'A1': A_l[0],
        'A2': A_l[1],
        'A3': A_l[2],
        'H0': getH0(lat, lon, time, nav_path)
        }
    return params

