# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 12:04:49 2021

@author: maria
"""
import numpy as np
import pandas as pd
import math
import datetime as datetime
from datetime import datetime as dtt

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

def GCtoLC(P0_cart, points_df):
    P0_g = cartToGeod(P0_cart[0], P0_cart[1], P0_cart[2])
    lat0 = degToRad(P0_g[0])
    lon0 = degToRad(P0_g[1])
    
    results_LC = pd.DataFrame()
    results_LC['dx'] = points_df['X'] - P0_cart[0]
    results_LC['dy'] = points_df['Y'] - P0_cart[1]
    results_LC['dz'] = points_df['Z'] - P0_cart[2]
    
    E=[]
    N=[]
    U=[]
    az = []
    el = []
    dist = []
    
    # Definition of the rotation matrix
    R = np.array([ [-np.sin(lon0), np.cos(lon0), 0],
                    [-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0), np.cos(lat0)],
                    [np.cos(lat0)*np.cos(lon0), np.cos(lat0)*np.sin(lon0), np.sin(lat0)]
                    ])
    
    for i in range(len(results_LC)):
        delta_array_i = np.array([[results_LC['dx'][i]], [results_LC['dy'][i]], [results_LC['dz'][i]]])
        ENU =np.dot(R, delta_array_i)
        E.append(ENU[0][0])
        N.append(ENU[1][0])
        U.append(ENU[2][0])
        d_h = np.sqrt(ENU[0][0]**2 + ENU[1][0]**2)
        dist.append(d_h)
        if d_h < 0.1:
            az.append(0)
            el.append(90)
        else:
            az.append(radToDeg(np.arctan(ENU[0][0]/ENU[1][0])))
            el.append(radToDeg(np.arctan2(ENU[2][0], d_h)))
    
    if 'datetime' in points_df.columns:
        results_LC['time'] = points_df['datetime']
    
    if 'time'in points_df.columns:
        results_LC['time'] = points_df['time']
        
    if 'sv' in points_df.columns:
        results_LC['sv'] = points_df['sv']
    
    results_LC['E'] = E
    results_LC['N'] = N
    results_LC['U'] = U
    results_LC['Az'] = az
    results_LC['El'] = el
    results_LC['Dist'] = dist
    
    results_LC = results_LC.drop(columns=['dx', 'dy', 'dz'])
    
    return results_LC

def fixDateTime(df):
    if 'time' in df.columns:
        conv_t = []
        for i in range(len(df)):
            t_i = df['time'][i]
            t_i_conv = dtt.strptime(t_i, '%Y-%m-%d %H:%M:%S')
            conv_t.append(t_i_conv)
        df = df.drop(columns=['time'])
        df['time'] = conv_t
        
        return df
    else:
        raise Exception('No time columns to convert in the DataFrame')