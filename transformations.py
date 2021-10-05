# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 12:04:49 2021

@author: maria
"""
import numpy as np
import math

def cartToGeod(x, y, z):
    a = 6378137
    e = 0.0818191908426215
    e2 = e*e
    b = a*(np.sqrt(1 - e2))
    eb2 = (a*a - b*b)/(b*b)
    # radius computation 
    r = np.sqrt(x*x + y*y)
    # longitude
    lon = np.arctan2(y, x)
    #latitute
    psi = np.arctan2(z, (r*np.sqrt(1-e2)))
    lat = np.arctan2((z+eb2*b*np.power(np.sin(psi), 3)), (r - e2*a*np.power(np.cos(psi), 3)))
    N = a/np.sqrt(1 - e2*np.power(np.sin(lat), 2))
    h = r/(np.cos(lat)) - N
    lon = radToDeg(lon)
    lat = radToDeg(lat)    
    return [lat, lon, h]

def geodToCart(lat, lon, h):
    a = 6378137
    e = 0.0818191908426215
    e2 = e*e
    lat = degToRad(lat)
    lon = degToRad(lon)
    N = a/np.sqrt(1 - e2*np.power(np.sin(lat), 2))
    x = (N + h)*np.cos(lat)*np.cos(lon)
    y = (N + h)*np.cos(lat)*np.sin(lon)
    z = (N*(1-e2) + h)*np.sin(lat)
    
    return [x, y, z]

def radToDeg(rad):
    deg = rad*180/math.pi
    return deg

def degToRad(deg):
    rad = deg*math.pi/180
    return rad

def topocent(XR, XS):
    XRgeod = cartToGeod(XR[0][0], XR[1][0], XR[2][0])
    phi = degToRad(XRgeod[0])
    lam = degToRad(XRgeod[1])
    
    cl = np.cos(lam)
    sl = np.sin(lam)
    cb = np.cos(phi)
    sb = np.sin(phi)
    F = np.array([[-sl, -sb*cl, cb*cl], [cl, -sb*sl, cb*sl], [0, cb, sb]])
    F_t = F.transpose()
    delta_X = XS - XR
    LC = np.dot(F_t, delta_X)
    E = LC[0][0]
    N = LC[1][0]
    U = LC[2][0]
    hor_dis = np.sqrt(E*E + N*N)
    if hor_dis < 0.1:
        Az = 0
        El = 90
    else:
        Az = radToDeg(np.arctan(E/N))
        El = radToDeg(np.arctan2(U,hor_dis))
    return [hor_dis, Az, El]

def topocent2(XR, XS):
    XScart = geodToCart(XS[0][0], XS[1][0], XS[2][0])
    XRgeod = cartToGeod(XR[0][0], XR[1][0], XR[2][0])
    phi = degToRad(XRgeod[0])
    lam = degToRad(XRgeod[1])
    XS_array = np.array([[XScart[0]], [XScart[1]], [XScart[2]]])
    
    cl = np.cos(lam)
    sl = np.sin(lam)
    cb = np.cos(phi)
    sb = np.sin(phi)
    F = np.array([[-sl, -sb*cl, cb*cl], [cl, -sb*sl, cb*sl], [0, cb, sb]])
    F_t = F.transpose()
    delta_X = XS_array - XR
    LC = np.dot(F_t, delta_X)
    E = LC[0][0]
    N = LC[1][0]
    U = LC[2][0]
    hor_dis = np.sqrt(E*E + N*N)
    if hor_dis < 0.1:
        Az = 0
        El = 90
    else:
        Az = (np.arctan(E/N))
        El = radToDeg(np.arctan2(U,hor_dis))
    return [hor_dis, Az, El]

def GCtoLC(P0, Px):
    P0_deg = cartToGeod(P0[0], P0[1], P0[2])
    lat0 = degToRad(P0_deg[0])
    lon0 = degToRad(P0_deg[1])
    delta_x = np.array([ [Px[0] - P0[0]],
                         [Px[1] - P0[1]],
                         [Px[2] - P0[2]] ])
    R = np.array([ [-np.sin(lon0), np.cos(lon0), 0],
                    [-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0), np.cos(lat0)],
                    [np.cos(lat0)*np.cos(lon0), np.cos(lat0)*np.sin(lon0), np.sin(lat0)]
                    ])
    LC = np.dot(R, delta_x)
    E = LC[0][0]
    N = LC[1][0]
    U = LC[2][0]
    return [E, N, U]

