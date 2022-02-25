# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 15:49:53 2022

@author: maria
"""

import functions as fn
import pandas as pd
import datetime as datetime
from datetime import datetime as dtt
import georinex as gr
import numpy as np
import transformations as trf 
import pointPositioning as pp
import easyPlot
import RINEXreader

import matplotlib.pyplot as plt
import datetime as datetime
from datetime import datetime as dtt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import matplotlib.patches as mpatches

sat = pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/sat_galileo_new.csv')

MOSE_cart = [4642432.701789316, 1028629.1051167124, 4236854.058403561]
MOSE_geod = trf.cartToGeod(4642432.701789316, 1028629.1051167124, 4236854.058403561)

elevation = []
azimuth = []
conv_t = []
for i in range(len(sat)):
    # Conversione tempo
    t_i = sat['time'][i]
    t_i_conv = dtt.strptime(t_i, '%Y-%m-%d %H:%M:%S')
    conv_t.append(t_i_conv)

sat_LC = trf.GCtoLC(MOSE_cart, sat)
for i in range(len(sat_LC)):
    E, N, U = sat_LC['E'][i], sat_LC['N'][i], sat_LC['U'][i]
    hor_dis = np.sqrt(E**2 + N**2)
    if hor_dis < 0.1:
        Az = 0
        El = 90
    else:
        Az = trf.radToDeg(np.arctan(E/N))
        El = trf.radToDeg(np.arctan2(U,hor_dis))
        
    azimuth.append(Az)
    elevation.append(El)

sat['El'] = elevation
sat['Az'] = azimuth    

sat = sat.drop(columns=['Unnamed: 0', 'time'])
sat['time'] = conv_t
sat['iono_delay'] = sat['TEC_new']*40.3*10**16/(1575.42*10**6)**2

nav_path = 'C:/git/pointPositioning/pointPositioning/mose_nav.rnx'
iono_params = RINEXreader.getIonoParams(nav_path, 'G')

klo = []
for i in range(len(sat)):
    w = sat['time'][i]
    start_day = w.day - w.weekday() - 1
    t = (w.day - start_day)*24*3600 + w.hour*3600 + w.minute*60 + w.second
    iono_i = pp.ionoCorrectionGPS(trf.radToDeg(MOSE_geod[0]), trf.radToDeg(MOSE_geod[1]), sat['Az'][i], sat['El'][i], t, iono_params)
    klo.append(iono_i)

sat['klobuchar'] = klo

sat_i = sat[sat['sv'] == 'E03'].reset_index()

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(sat_i['El'], sat_i['iono_delay'], 'o', color='blue')
ax.plot(sat_i['El'], sat_i['klobuchar'], 'o', color='red')
ax.set(xlabel='Elevation', ylabel='Ionospheric delay', title = sat_i['sv'][0])
o_patch = mpatches.Patch(color='blue', label='NeQuikG')
p_patch = mpatches.Patch(color='red', label='Klobuchar')
#ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
#ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
plt.legend(handles=[p_patch, o_patch])
plt.show()

