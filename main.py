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
import easyPlot

nav_path = 'C:/git/pointPositioning/pointPositioning/mose_nav.rnx'
obs_path = 'C:/git/pointPositioning/pointPositioning/mose_obs.crx'
eph_path = 'C:/git/pointPositioning/pointPositioning/EPH_2021_5_6.SP3'

time_range = []
START = dtt.strptime('2021-05-06 00:00:00', '%Y-%m-%d %H:%M:%S')
END = dtt.strptime('2021-05-06 23:45:00', '%Y-%m-%d %H:%M:%S')
t = START
while t <= END:
    time_range.append(t)
    t = t + datetime.timedelta(minutes=5)

# Calcolo delle orbite
'''
sat_galileo = fn.getGALILEOorbits(nav_path, obs_path, time_range)
sat_galileo.to_csv('C:/git/pointPositioning/pointPositioning/csv/sat_galileo.csv')
sat_gps = fn.getGPSorbits(nav_path, obs_path, time_range)
sat_gps.to_csv('C:/git/pointPositioning/pointPositioning/csv/sat_gps.csv')
satellites = pd.concat([sat_gps, sat_galileo])
satellites.to_csv('C:/git/pointPositioning/pointPositioning/csv/satellites.csv')
'''
sat_galileo = pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/sat_galileo.csv')
sat_galileo_new = pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/sat_galileo_new.csv')
sat_gps = pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/sat_gps.csv')
satellites = pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/satellites.csv')

conv_t = []
for i in range(len(sat_gps)):
    t_i = sat_gps['time'][i]
    t_i_conv = dtt.strptime(t_i, '%Y-%m-%d %H:%M:%S')
    conv_t.append(t_i_conv)

sat_gps = sat_gps.drop(columns=['Unnamed: 0', 'time'])
sat_gps['time'] = conv_t
sat_galileo_new['iono_delay'] = sat_galileo_new['TEC_new']*40.3*10**16/(1575.42*10**6)**2

klo = []

for i in range(len(sat_galileo_new)):
    

satellites = pd.concat([sat_gps, sat_galileo_new])
satellites = satellites.reset_index().drop(columns=['index'])
satellites['iono_delay'] = satellites['TEC']*40.3*10**16/(1575.42*10**6)**2

sat_galileo_new = sat_galileo_new.drop(columns=['time', 'Unnamed: 0'])
sat_galileo_new['time'] = conv_t
sat_galileo_new['iono_delay'] = sat_galileo_new['TEC']*40.3*10**16/(1575.42*10**6)**2

# Controllo delle posizioni e clock offset
'''
check_pos_gal = fn.checkSatPos(sat_galileo, eph_path)
check_pos_gps = fn.checkSatPos(sat_gps, eph_path)
'''

# Controllo delle velocità
'''
check_vel_GAL = fn.checkSatVelGAL(sat_galileo, nav_path, time_range)
check_vel_GPS = fn.checkSatVelGPS(sat_gps, nav_path, time_range)
'''

# POINT POSITIONING

cutoff = 5
results = pp.pointPositioning2(satellites, nav_path, obs_path, cutoff)
results_k = pp.pointPositioning2(satellites, nav_path, obs_path, cutoff)


results_gps = pp.pointPositioning(sat_gps, nav_path, obs_path, cutoff)

results_galileo = pp.pointPositioning(sat_galileo, nav_path, obs_path, cutoff)
results_galileo_with_iono = pp.pointPositioning3(sat_galileo_new, nav_path, obs_path, cutoff)



M0SE_cart = [4642432.701789316, 1028629.1051167124, 4236854.058403561]

results_LC = trf.GCtoLC(M0SE_cart, results)
results_LC_k = trf.GCtoLC(M0SE_cart, results_k)
results_LC_GPS = trf.GCtoLC(M0SE_cart, results_gps)
results_LC_GAL = trf.GCtoLC(M0SE_cart, results_galileo)
results_LC_GAL_iono = trf.GCtoLC(M0SE_cart, results_galileo_with_iono)

easyPlot.getPlot(results_LC, 'time', 'E', 'red')
easyPlot.getPlot(results_LC, 'time', 'N', 'blue')
easyPlot.getPlot(results_LC, 'time', 'U', 'green')
easyPlot.getPlot(results, 'datetime', 'dtr_GPS', 'magenta')

easyPlot.getDoublePlot(results, results_galileo, 'datetime', 'dtr')

easyPlot.getDoublePlot(results_LC_k, results_LC, 'time', 'E')
easyPlot.getDoublePlot(results_LC_k, results_LC, 'time', 'N')
easyPlot.getDoublePlot(results_LC_k, results_LC, 'time', 'U')

easyPlot.getTriplePlot(results_LC, results_LC_GPS, results_LC_GAL, 'time', 'E')
easyPlot.getTriplePlot(results_LC, results_LC_GPS, results_LC_GAL, 'time', 'N')
easyPlot.getTriplePlot(results_LC, results_LC_GPS, results_LC_GAL, 'time', 'U')

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(results['datetime'], results['dtr_GPS'], '-', color='orange')
ax.plot(results['datetime'], results['dtr_GAL'], '-', color='purple')
#ax.plot(results['datetime'], results['dtr_GAL']-results['dtr_GPS'], '-', color='red')
ax.set(xlabel='time', ylabel='dtr', title = 'Time-dtr')
o_patch = mpatches.Patch(color='orange', label='dtr GPS')
p_patch = mpatches.Patch(color='purple', label='dtr GAL')
ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
plt.legend(handles=[p_patch, o_patch])
plt.show()




# Analsis on crytical epochs

t1 = results_LC_GAL[results_LC_GAL['E'] == abs(results_LC_GAL['E']).max()].reset_index()['time'][0]

t2 = results_LC_GAL[results_LC_GAL['N'] == abs(results_LC_GAL['N']).max()].reset_index()['time'][0]

t3 = results_LC_GAL[results_LC_GAL['U'] == abs(results_LC_GAL['U']).max()].reset_index()['time'][0]

mask_sat = sat_galileo[sat_galileo['time'] == t1]

'''TEST TIME SYSTEM CORR'''
GAGP = [3.2014213502*10**-9, -2.220446049*10**-15, 345600, 2156]

time_deltas = []

for i in time_range:
    start_day = i.day - i.weekday() - 1
    t = (i.day - start_day)*24*3600 + i.hour*3600 + i.minute*60 + i.second
    delta_t_i = GAGP[0] + GAGP[1]*(t - GAGP[2])
    time_deltas.append(delta_t_i)

time_deltas_s = pd.DataFrame(time_deltas, columns=['GAGP'])

results_galileo['dtr_GAGP'] = results_gps['dtr'] + time_deltas_s['GAGP']

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(results_galileo['datetime'], (results_galileo['dtr']-results_gps['dtr']), '-', color='orange')
#ax.plot(results_galileo['datetime'], results_galileo['dtr'], '-', color='purple')
ax.set(xlabel='time', ylabel='dtr', title = '')
#o_patch = mpatches.Patch(color='orange', label='dtr_GPS + GAGP')
#p_patch = mpatches.Patch(color='purple', label='GAL')
ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
#plt.legend(handles=[p_patch, o_patch])
plt.show()

