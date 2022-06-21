# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 10:16:39 2022

@author: maria
"""
# Import of the required libraries
import pandas as pd
import datetime as datetime
from datetime import datetime as dtt
import numpy as np

# Import of the software's modules
import codepos.RINEXreader as rr
import codepos.easyPlot as ep
import codepos.functions as fn
import codepos.pointPositioning as pp
import codepos.transformations as trf
from NQ import TEC as modTEC
from NQ import modelParams as mp

'''1) path to RINEX files'''
day = '14'
path = 'C:/git/pointPositioning/pointPositioning/testRINEX/2021_07_13/'
#path = 'C:/git/pointPositioning/pointPositioning/MN/'+day+'-07-2021/'
nav_path = path + 'mil.rnx'
nav_path2 = path + 'nav2.rnx'
nav_path3 = path + 'nav3.rnx'
nav_path_brdc = path + 'nav_brdc.rnx'
nav_path_mos = path + 'nav_m0se.rnx'
obs_path = path + 'obs.crx'
eph_path = path + 'eph.SP3'

'''2) Definition of the time range'''
time_range = []
START = dtt.strptime('2021-07-'+day+' 00:00:00', '%Y-%m-%d %H:%M:%S')
END = dtt.strptime('2021-07-'+day+' 23:45:00', '%Y-%m-%d %H:%M:%S')
t = START
while t <= END:
    time_range.append(t)
    t = t + datetime.timedelta(minutes=5)

'''3) Coordinates of the receiver'''
M0SE_cart = [4642432.701789316, 1028629.1051167124, 4236854.058403561]
M0SE_geod = trf.cartToGeod(4642432.701789316, 1028629.1051167124, 4236854.058403561)

MIL_geod = [45.478368595, 9.229212847777777, 191.125]
MIL_cart = trf.geodToCart(45.478368595, 9.229212847777777, 191.125)

'''4) Orbit computation'''

satellites = fn.getBRDCorbits(nav_path_brdc, time_range, ['G', 'E'])

# sat_galileo_mos = fn.getGALILEOorbits(nav_path_mos, obs_path, time_range) #56 secondi
# sat_gps_mos = fn.getOrbits(nav_path_mos, obs_path, time_range)

# sat_galileo_brdc = fn.getGALILEOorbits(nav_path_brdc, obs_path, time_range) #56 secondi
# sat_gps_brdc = fn.getOrbits(nav_path_brdc, obs_path, time_range)

# sat_galileo = fn.getGALILEOorbits(nav_path, obs_path, time_range) #56 secondi
# sat_gps = fn.getGPSorbits(nav_path, obs_path, time_range)

# sat_galileo2 = fn.getGALILEOorbits(nav_path2, obs_path, time_range) #56 secondi
# sat_gps2 = fn.getGPSorbits(nav_path2, obs_path, time_range)

# sat_galileo3 = fn.getGALILEOorbits(nav_path3, obs_path, time_range) #56 secondi
# sat_gps3 = fn.getGPSorbits(nav_path3, obs_path, time_range)

# sat_galileo = trf.fixDateTime(pd.read_csv(path+'sat_gal.csv'))
# sat_gps = trf.fixDateTime(pd.read_csv(path+'sat_gps.csv'))
# satellites = trf.fixDateTime(pd.read_csv(path+'satellites.csv'))

'''5) Check of satellites positions and clock offsets''' 

check_pos = fn.checkSatPos(satellites, eph_path)

# check_pos_gal_brdc = fn.checkSatPos(sat_galileo_brdc, eph_path)
# check_pos_gps_brdc = fn.checkSatPos(sat_gps_brdc, eph_path)

# check_temp = check_pos_gps.query('delta_x < 5')
# check_temp = check_temp.query('delta_y < 5')
# check_temp = check_temp.query('delta_z < 5')

# av_gps_sat = sat_gps.iloc[check_temp.index]
# sat_gps = av_gps_sat.reset_index().drop(columns=['index'])

'''6) Galileo Ionospheric delay with NeQuickG algorithm'''

parametri = []
for time in time_range:
    p = mp.getModelParams(M0SE_geod[0], M0SE_geod[1], time, nav_path_brdc)
    parametri.append(p)

dict_par = dict(zip(time_range, parametri))
    
P1 = M0SE_geod
P1[2] = P1[2]*10**(-3)


delay = []
for i in range(len(sat_galileo)):
    time = sat_galileo['time'][i]
    par = dict_par[time]
    PS = trf.cartToGeod(sat_galileo['xs'][i], sat_galileo['ys'][i], sat_galileo['zs'][i])
    PS[2] = PS[2]*10**(-3)
    d = modTEC.getTEC(P1, PS, par)
    delay.append(d*40.3*10**16/(1575.42*10**6)**2)

sat_galileo_brdc['iono_delay'] = delay
sat_gps_brdc['iono_delay'] = np.nan

'''7) Merge of the satellites dataframe'''
satellites = pd.concat([sat_gps_brdc, sat_galileo_brdc])
satellites.to_csv(path + 'satellites.csv')

'''8) Check on velocities -> TO BE DEFINED'''
# check_vel_GAL = fn.checkSatVelGAL(sat_galileo, nav_path, time_range)
# check_vel_GPS = fn.checkSatVelGPS(sat_gps, nav_path, time_range)

'''9) POINT POSITIONING'''

cutoff = 5
results = pp.pointPositioning(satellites, nav_path_brdc, obs_path, cutoff)

# Difference of the estimated clock offsets
results['dtr'] = results['dtr_GAL'] - results['dtr_GPS']

results.to_csv(path + 'results.csv')

'''10) Conversion in Local Cartesian'''
results = results.rename(columns={'xr':'X', 'yr':'Y', 'zr':'Z'})
results_LC = trf.GCtoLC(M0SE_cart, results)

results_LC.to_csv(path + 'results_LC.csv')


'''11) Plotting'''
ep.getPlot(results_LC, 'time', 'E', 'red')
ep.getPlot(results_LC, 'time', 'N', 'blue')
ep.getPlot(results_LC, 'time', 'U', 'green')
ep.getPlot(results, 'datetime', 'dtr', 'orange')


import matplotlib.pyplot as plt
import datetime as datetime
from datetime import datetime as dtt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import matplotlib.patches as mpatches

fig, ax = plt.subplots(figsize=(10,6))
ax.plot(results['datetime'], results['dtr_GAL'], '-', color='purple')
ax.plot(results['datetime'], results['dtr_GPS'], '-', color='pink')
ax.set(xlabel='time', ylabel='dtr', title = 'Time - dtr')
o_patch = mpatches.Patch(color='purple', label='GAL')
p_patch = mpatches.Patch(color='pink', label='GPS')
ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
plt.legend(handles=[p_patch, o_patch])
plt.show()

deltaDTR = pd.DataFrame()
delta_list = []
for i in range(len(results)-1):
    Ddtr = results['dtr_GAL'][i+1] - results['dtr_GAL'][i]
    delta_list.append(Ddtr)

deltaDTR['Ddtr_gal'] = delta_list
deltaDTR['time'] = time_range[0:(len(time_range)-1)]

ep.getPlot(deltaDTR, 'time', 'Ddtr_gal', 'red')

results_k = pp.pointPositioning2(satellites, nav_path, obs_path, cutoff)

results_gal_NQ = pp.pointPositioning3(sat_galileo, nav_path, obs_path, cutoff)
results_gal_k = pp.pointPositioning3(sat_galileo, nav_path, obs_path, cutoff)
results_gal_0 = pp.pointPositioning3(sat_galileo, nav_path, obs_path, cutoff)

results.to_csv('C:/git/pointPositioning/pointPositioning/csv/pos/results.csv')
results_k.to_csv('C:/git/pointPositioning/pointPositioning/csv/pos/results_k.csv')
results_gal_NQ.to_csv('C:/git/pointPositioning/pointPositioning/csv/pos/results_gal_NQ.csv')
results_gal_k.to_csv('C:/git/pointPositioning/pointPositioning/csv/pos/results_gal_k.csv')
results_gal_0.to_csv('C:/git/pointPositioning/pointPositioning/csv/pos/results_gal_0.csv')

M0SE_cart = [4642432.701789316, 1028629.1051167124, 4236854.058403561]

results_LC = trf.GCtoLC(M0SE_cart, results)
results_LC_k = trf.GCtoLC(M0SE_cart, results_k)
results_LC_gal = trf.GCtoLC(M0SE_cart, results_gal_NQ)
results_LC_gal_k = trf.GCtoLC(M0SE_cart, results_gal_k)
results_LC_gal_0 = trf.GCtoLC(M0SE_cart, results_gal_0)

ep.getPlot(results_LC, 'time', 'E', 'red')
ep.getPlot(results_LC, 'time', 'N', 'blue')
ep.getPlot(results_LC, 'time', 'U', 'green')
easyPlot.getPlot(results, 'datetime', 'dtr_GPS', 'magenta')

# easyPlot.getDoublePlot(results, results_galileo, 'datetime', 'dtr')

# easyPlot.getDoublePlot(results_LC_k, results_LC, 'time', 'E')
# easyPlot.getDoublePlot(results_LC_k, results_LC, 'time', 'N')
# easyPlot.getDoublePlot(results_LC_k, results_LC, 'time', 'U')

# easyPlot.getTriplePlot(results_LC_gal, results_LC_gal_k, results_LC_gal_0, 'time', 'E')
# easyPlot.getTriplePlot(results_LC_gal, results_LC_gal_k, results_LC_gal_0, 'time', 'N')
# easyPlot.getTriplePlot(results_LC_gal, results_LC_gal_k, results_LC_gal_0, 'time', 'U')

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

