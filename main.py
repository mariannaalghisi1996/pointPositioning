# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 10:16:39 2022

@author: maria
"""
# Import of the required libraries
import pandas as pd
import datetime as datetime
from datetime import datetime as dtt

# Import of the software's modules
import codepos.RINEXreader as rr
import codepos.easyPlot as ep
import codepos.functions as fn
import codepos.pointPositioning as pp
import codepos.transformations as trf
from NQ import TEC

nav_path = 'C:/git/pointPositioning/pointPositioning/testRINEX/mose_nav.rnx'
obs_path = 'C:/git/pointPositioning/pointPositioning/testRINEX/mose_obs.crx'
eph_path = 'C:/git/pointPositioning/pointPositioning/testRINEX/EPH_2021_5_6.SP3'

time_range = []
START = dtt.strptime('2021-05-06 00:00:00', '%Y-%m-%d %H:%M:%S')
END = dtt.strptime('2021-05-06 23:55:00', '%Y-%m-%d %H:%M:%S')
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
sat_galileo = trf.fixDateTime(pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/orbite/sat_galileo.csv'))
sat_gps = trf.fixDateTime(pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/orbite/sat_gps.csv'))
satellites = trf.fixDateTime(pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/orbite/satellites.csv'))

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

# Plot delle orbite
import folium

M0SE_cart = [4642432.701789316, 1028629.1051167124, 4236854.058403561]
M0SE_geod = trf.cartToGeod(4642432.701789316, 1028629.1051167124, 4236854.058403561)

m = folium.Map(location=[M0SE_geod[0], M0SE_geod[1]], zoom_start=5)

G1 = sat_gps[sat_gps['sv'] == 'G01'].reset_index().drop(columns=['index'])
lat, lon = [], []
for i in range(len(G1)):
    temp = trf.cartToGeod(G1['xs'][i], G1['ys'][i], G1['zs'][i])
    folium.Circle(
        location=(temp[0], temp[1]),
        radius=10,
        fill=True,
        fill_opacity=0.7
    ).add_to(m)

m.save('map.html')


# POINT POSITIONING

cutoff = 5
results = pp.pointPositioning2(satellites, nav_path, obs_path, cutoff)
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

easyPlot.getPlot(results_LC_gal_k, 'time', 'E', 'red')
easyPlot.getPlot(results_LC_k, 'time', 'N', 'blue')
easyPlot.getPlot(results_LC_k, 'time', 'U', 'green')
easyPlot.getPlot(results, 'datetime', 'dtr_GPS', 'magenta')

easyPlot.getDoublePlot(results, results_galileo, 'datetime', 'dtr')

easyPlot.getDoublePlot(results_LC_k, results_LC, 'time', 'E')
easyPlot.getDoublePlot(results_LC_k, results_LC, 'time', 'N')
easyPlot.getDoublePlot(results_LC_k, results_LC, 'time', 'U')

easyPlot.getTriplePlot(results_LC_gal, results_LC_gal_k, results_LC_gal_0, 'time', 'E')
easyPlot.getTriplePlot(results_LC_gal, results_LC_gal_k, results_LC_gal_0, 'time', 'N')
easyPlot.getTriplePlot(results_LC_gal, results_LC_gal_k, results_LC_gal_0, 'time', 'U')

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

