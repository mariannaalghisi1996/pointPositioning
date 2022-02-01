# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 10:16:39 2022

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

nav_path = 'C:/git/pointPositioning/pointPositioning/mose_nav.rnx'
obs_path = 'C:/git/pointPositioning/pointPositioning/mose_obs.crx'

time_range = []
START = dtt.strptime('2021-05-06 00:00:00', '%Y-%m-%d %H:%M:%S')
END = dtt.strptime('2021-05-06 23:45:00', '%Y-%m-%d %H:%M:%S')
t = START
while t <= END:
    time_range.append(t)
    t = t + datetime.timedelta(minutes=5)

sat_galileo = fn.getGALILEOorbits(nav_path, obs_path, time_range)

sat_gps = fn.getGPSorbits(nav_path, obs_path, time_range)

satellites = pd.concat([sat_gps, sat_galileo])

cutoff = 5
results = pp.pointPositioning(satellites, nav_path, obs_path, cutoff)

results_gps = pp.pointPositioning(sat_gps, nav_path, obs_path, cutoff)
results_galileo = pp.pointPositioning(sat_galileo, nav_path, obs_path, cutoff)



M0SE_cart = [4642432.701789316, 1028629.1051167124, 4236854.058403561]

results_LC = trf.GCtoLC(M0SE_cart, results)
results_LC_GPS = trf.GCtoLC(M0SE_cart, results_gps)
results_LC_GAL = trf.GCtoLC(M0SE_cart, results_galileo)

import matplotlib.pyplot as plt
import datetime as datetime
from datetime import datetime as dtt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import matplotlib.patches as mpatches

print('MEDIA ERRORI:')
print('Est: ', results_LC['E'].mean(), '\n',
      'Nord: ', results_LC['N'].mean(), '\n',
      'Quota: ', results_LC['U'].mean() )

print('DEVIAZIONE STANDARD:', '\n',
      'Est: ', results_LC['E'].std(), '\n',
      'Nord: ', results_LC['N'].std(), '\n',
      'Quota: ', results_LC['U'].std() )

print('Valori Massimi:', '\n',
      'Est: ', abs(results_LC['E']).max(), '\n',
      'Nord: ', abs(results_LC['N']).max(), '\n',
      'Quota: ', abs(results_LC['U']).max() )

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(results_LC['time'], results_LC['U'],
        '-',
        color='orange')
ax.plot(results_LC_no_el['time'], results_LC_no_el['U'],
        '-',
        color='purple')
ax.set(xlabel="Time", ylabel='UP', title = 'Time - Differences in UP')
ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
o_patch = mpatches.Patch(color='orange', label='GPS+GAL')
p_patch = mpatches.Patch(color='purple', label='GPS+GAL no elevazione')
plt.legend(handles=[o_patch, p_patch])
plt.show()

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(results_LC_no_el['time'], results_LC['N'],
        '-',
        color='blue')
ax.set(xlabel="Time", ylabel='NORTH', title = 'Time-NORTH')
ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
plt.show()

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(results_LC['time'], results_LC['U'],
        '-',
        color='green')
ax.set(xlabel="Time", ylabel='UP', title = 'Time-UP')
ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
plt.show()


fig, ax = plt.subplots(figsize=(10,6))
ax.plot(results['datetime'], abs(results_galileo['dtr'] - results_gps['dtr']),
        '-',
        color='purple')
ax.set(xlabel="Time", ylabel='dtr', title = 'Time-dtr')
ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
plt.show()



test = results_LC[results_LC['U'] == results_LC['U'].max()]
time_max = dtt.strptime('2021-05-06 09:05:00', '%Y-%m-%d %H:%M:%S')

sat_time_max = satellites[satellites['time']==time_max].reset_index()
