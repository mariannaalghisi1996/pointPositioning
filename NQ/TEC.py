# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 15:25:56 2022

@author: maria
"""
import numpy as np
import math
import georinex as gr

import NQ.electronDensity as ed
import NQ.modelParams as mp
from NQ.constants import R_E, DR, RD


def vertGNtest(h1, h2, n, p):
    delta_n = (h2 - h1)/n
    g = 0.5773502691896*delta_n
    y = g + (delta_n - g)/2
    GN2 = 0
    for i in range(n):
        hb = y + i*delta_n
        ht = y + i*delta_n + g
        GN2 = GN2 + ed.getN(hb, p) + ed.getN(ht, p)
    GN2 = GN2*delta_n/2
    
    GN1 = GN2
    n = 2*n
    delta_n = (h2 - h1)/n
    g = 0.5773502691896*delta_n
    y = g + (delta_n - g)/2
    GN2 = 0
    for i in range(n):
        hb = y + i*delta_n
        ht = y + i*delta_n + g
        GN2 = GN2 + ed.getN(hb, p) + ed.getN(ht, p)
    GN2 = GN2*delta_n/2
    
    return [GN1, GN2]

def verticalTEC(h1, h2, eps, p):
    n = 8
    GN = vertGNtest(h1, h2, n, p)
    GN1 = GN[0]
    GN2 = GN[1]
    
    while(abs(GN2 - GN1) > eps*abs(GN1)):
        n = 2*n
        GN = vertGNtest(h1, h2, n, p)
        GN1 = GN[0]
        GN2 = GN[1]
    
    vTEC = (GN2 + (GN2 - GN1)/15)*10**(-13)
    
    return vTEC

def getZenit(lat1, lon1, h1, lat2, lon2, h2):
    '''
    Returns Zenith angle between P1 and P2
    Latitude and longitude values in iputs are in radians
    Heights are in kilometers
    '''
    r1 = h1 + R_E
    r2 = h2 + R_E
    cosZ = np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(lon2 - lon1)
    sinZ = np.sqrt(1 - cosZ**2)
    Z = np.arctan2(sinZ, cosZ - r1/r2)
    
    return [Z, cosZ, sinZ]
  

def slantGeometricalConf(lat1, lon1, h1, lat2, lon2, h2):
    lat1 = lat1*DR
    lon1 = lon1*DR
    lat2 = lat2*DR
    lon2 = lon2*DR
    r1 = h1 + R_E
    r2 = h2 + R_E
    
    # Zenith angle computation:
    # Z is the Earth angle on the great circle connecting the receiver P1 and the satellite P2
    cosD = np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(lon2 - lon1)
    sinD = np.sqrt(1 - cosD**2)
    Z = np.arctan2(sinD, cosD - r1/r2)
    
    # Ray-perigee computation Pp:
    # Ray-perigee radius
    rp = r1*np.sin(Z)
    # Ray-perigee latitude latp:
    if abs(abs(lat1*RD) - 90) < 10**(-10):
        if lat1*RD > 0:
            latp = Z
        else:
            latp = -Z
    else:
        sinS = (np.sin(lon2 - lon1)*np.cos(lat2))/sinD
        cosS = (np.sin(lat2) - cosD*np.sin(lat1))/(sinD*np.cos(lat1))
        Dp = math.pi/2 - Z
        sin_latp = np.sin(lat1)*np.cos(Dp) - np.cos(lat1)*np.sin(Dp)*cosS
        cos_latp = np.sqrt(1 - sin_latp**2)
        latp = np.arctan2(sin_latp, cos_latp)
    # Ray-perigee longitude lonp:
    if abs(abs(lat1*RD) - 90) < 10**(-10):
        if lat1*RD >= 0:
            lonp = lon2 + math.pi
        else:
            lonp = lon2
    else:
        sin_lonplon1 =-(sinS*np.sin(Dp))/cos_latp
        cos_lonplon1 = (np.cos(Dp) - np.sin(lat1)*sin_latp)/(np.cos(lat1)*cos_latp)
        lonp = np.arctan2(sin_lonplon1, cos_lonplon1) + lon1
    
    # Great circle angle GCA from ray-perigee to satellite:
    if abs(abs(latp*RD)-90) < 10**(-10):
        W = abs(lat2 - latp)
    else:
        cosW = sin_latp*np.sin(lat2) + cos_latp*np.cos(lat2)*np.cos(lon2 - lonp)
        sinW = np.sqrt(1 - cosW**2)
        W = np.arctan2(sinW, cosW)
    
    # Sin and cos of azimuth of sat as seen from ray-perigee Pp:
    if abs(abs(latp*RD)-90) < 10**(-10):
        sinSp = 0
        if latp > 0:
            cosSp = -1
        else:
            cosSp = 1
    else:
        sinSp = np.cos(lat2)*np.sin(lon2 - latp)/sinW
        cosSp = (np.sin(lat2) - sin_latp*cosW)/(cos_latp*sinW)
    
    geomDict = {
        'rp': rp,
        'latp': latp,
        'lonp': lonp,
        'sinSp': sinSp,
        'cosSp': cosSp,
        'W': W}   
    
    return geomDict

def c(s, rpParams):
    '''
    Function that returns the coordinates of a point with a given distance s 
    with respect to the perigee.

    Parameters
    ----------
    s : distance wrt perigee (km)
    rpParams : ray-perigee parameters

    Returns
    -------
    [lats, lons, hs] with lats, lons in radians and hs in kilometers
    '''
    rp, latp, lonp = rpParams.get('rp'), rpParams.get('latp'), rpParams.get('lonp')
    cosSp, sinSp = rpParams.get('cosSp'), rpParams.get('sinSp')
    
    hs = np.sqrt(s**2 + rp**2) - R_E
    
    # Computation of great circle parameters
    tanDs = s/rp
    cosDs = 1/np.sqrt(1 + tanDs**2)
    sinDs = tanDs*cosDs
    
    sin_lats = np.sin(latp)*cosDs + np.cos(latp)*sinDs*cosSp
    cos_lats = np.sqrt(1 - sin_lats**2)
    lats = np.arctan2(sin_lats, cos_lats)
    
    sin_lonslonp = sinDs*sinSp*np.cos(latp)
    cos_lonslonp = cosDs - np.sin(latp)*sin_lats
    lons = (np.arctan2(cos_lonslonp, sin_lonslonp) + lonp)
    
    return [lats, lons, hs]
    

def slantGNtest(g1, g2, n, p, rpParams):
    print('n: ', n)
    delta_n = (g2 - g1)/n
    g = 0.5773502691896*delta_n
    y = g1 + (delta_n - g)/2
    GN2 = 0
    for i in range(n):
        n_coord_1 = c(y + i*delta_n, rpParams)
        n_coord_2 = c(y + i*delta_n + g, rpParams)
        hs_1 = n_coord_1[2]
        hs_2 = n_coord_2[2]
        #print('Round ', i, ' h: ', hs_1, hs_2)
        Summ_i = ed.getN(hs_1, p) + ed.getN(hs_2, p)
        GN2 = GN2 + Summ_i
    GN2 = GN2*delta_n/2
    
    n = n*2
    GN1 = GN2
    GN2 = 0
    for i in range(n):
        n_coord_1 = c(y + i*delta_n, rpParams)
        n_coord_2 = c(y + i*delta_n + g, rpParams)
        hs_1 = n_coord_1[2]
        hs_2 = n_coord_2[2]
        Summ_i = ed.getN(hs_1, p) + ed.getN(hs_2, p)
        GN2 = GN2 + Summ_i
    GN2 = GN2*delta_n/2
    
    return [GN1, GN2]

def gauss(g1, g2, eps, rpParams, p):
    '''
    Parameters
    ----------
    g1 : distances from the ray perigee of the first integration endpoint [km].
    g2 : distances from the ray perigee of the second integration endpoint [km].
    eps : target integration accuracy (max: 10**(-3)).
    rpParams : rp, sin(latp), cos(latp), sinOp, cosOp, lonp.
    alpha : Azimuth parameters a0, a1, a2 transmitted by rinex navigation file
    mth : month
    UT : Universal time

    Returns
    -------
    None.
    '''
    n = 8
    
    GN = slantGNtest(g1, g2, n, p, rpParams)
    GN1 = GN[0]
    GN2 = GN[1]
    
    n = n*2
    GN_next = slantGNtest(g1, g2, n, p, rpParams)
    GN1_next = GN[0]
    GN2_next = GN[1]
    
    while(abs(GN1 - GN1_next)/abs(GN1_next) > eps) or (abs(GN2 - GN2_next)/abs(GN2_next) > eps):
        n = 2*n
        GN1 = GN1_next
        GN2 = GN2_next
        GN_next = slantGNtest(g1, g2, n, p, rpParams)
        GN1_next = GN_next[0]
        GN2_next = GN_next[1]
    
    TEC = (GN2 - (GN2 - GN1)/15)*10**(-13)
    
    return TEC   
        
    
    
def slantTEC(P1, P2, p, rpParams):
    '''
    P1, P2 = [lat, lon, h(km)]

    Parameters
    ----------
    P1 : coord of the receiver
    P2 : coord of the satellite
    time : Y/m/d HH:MM:SS
    nav_path : Path to navigation file 
    '''
    lat1, lon1, h1 = P1[0]*DR, P1[1]*DR, P1[2]
    lat2, lon2, h2 = P2[0]*DR, P2[1]*DR, P2[2]
    r1 = h1 + R_E
    r2 = h2 + R_E
    # Geometrical parameters
    rp = rpParams.get('rp')
    # Integration endpoints s1, s2 (distances of P1 and P2 respectively from the ray perigee):
    s1 = np.sqrt(r1**2 - rp**2)
    s2 = np.sqrt(r2**2 - rp**2)
    
    
    if (h1<1000 and h2>2000):
        sa = np.sqrt((R_E + 1000)**2 - rp**2)
        sb = np.sqrt((R_E + 2000)**2 - rp**2)
        sTEC = gauss(s1, sa, 0.001, rpParams, p) + gauss(sa, sb, 0.01, rpParams, p) + gauss(sb, s2, 0.01, rpParams, p)
    else:
        sTEC = gauss(s1, s2, 0.01, rpParams, p)
    
    return sTEC

def getTEC(P1, P2, params):
    lat1, lon1, h1 = P1[0], P1[1], P1[2]
    lat2, lon2, h2 = P2[0], P2[1], P2[2]
    rpParams = slantGeometricalConf(lat1, lon1, h1, lat2, lon2, h2)
    rp = rpParams.get('rp')
    
    if rp < 0.1:
        print('Vertical TEC algorithm used')
        if h2 > 1000:
            eps = 0.01
        else:
            eps = 0.001
        TEC = verticalTEC(h1, h2, eps, params)
    else:
        print('Slant TEC algorithm used')
        TEC = slantTEC(P1, P2, params, rpParams)
    
    return TEC
    
