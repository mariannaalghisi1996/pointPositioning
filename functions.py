# -*- coding: utf-8 -*-
"""
@author: Marianna
"""

import datetime as datetime
from datetime import datetime as dtt
import pandas as pd
import numpy as np
import RINEXreader
import georinex as gr

def fixTime(df):
    ID_list = []
    maxx = df['time'].max()
    for i in range(len(df)):
        t = df['time'][i]
        while t.second != 0:
            t = t + datetime.timedelta(seconds = 16)
        df['time'][i] = t
        if t.day != maxx.day:
            ID_list.append(t.hour)
        else:
            ID_list.append(24)
    df['ID'] = ID_list
    
    return df

def getPar(sv, time, df_gps):
    mask = df_gps['sv'] == sv
    df_gps_2 = df_gps[mask]
    id_list = df_gps_2['ID'].to_list()
    IDk = time.hour
    if (IDk in id_list):
        params_m = df_gps_2['ID'] == IDk
        params = (df_gps_2[params_m]).reset_index().drop(columns=['index', 'time'])
        if type(params) != None:
            params['time'] = time
            return params
    else:
        if ((IDk-1) in id_list):
            params_m = df_gps_2['ID'] == (IDk-1)
            params = (df_gps_2[params_m]).reset_index().drop(columns=['index', 'time'])
            if type(params) != None:
                params['time'] = time
                return params
        else: 
            if ((IDk+1) in id_list):
                params_m = df_gps_2['ID'] == (IDk+1)
                params = (df_gps_2[params_m]).reset_index().drop(columns=['index', 'time'])
                if type(params) != None:
                    params['time'] = time
                    return params
            else:
                if ((IDk+2) in id_list):
                    params_m = df_gps_2['ID'] == (IDk + 2)
                    params = (df_gps_2[params_m]).reset_index().drop(columns=['index', 'time'])
                    if type(params) != None:
                        params['time'] = time
                        return params
        return 'error'  
    
def satPosVel(res):
    mu =  398600500000000 #m3/sec2
    omega_dot_E = 7.2921151467*0.00001 #rad/sec
    i = 0
    time = res['time'][i]
    start_day = time.day - time.weekday() - 1
    t = (time.day - start_day)*24*3600 + time.hour*3600 + time.minute*60 + time.second
    # Compute the time tk from the ephemerides reference epoch toe 
    tk = t - res['Toe'][i]
    if tk> 302400:
        tk = tk - 604800
    elif tk <-302400:
         tk = tk + 604800
    dts = res['SVclockBias'][i] + tk*res['SVclockDrift'][i] + tk*tk*res['SVclockDriftRate'][i]
    
    """SATELLITE POSITION COMPUTATION"""
    # Semi-major axis
    A = res['sqrtA'][i]*res['sqrtA'][i]
    # Computed mean motion
    n0 = np.sqrt(mu/(np.power(A, 3)))
    # Corrected mean motion
    n = n0 + res['DeltaN'][i]
    # Compute the mean anomaly for tk
    Mk = res['M0'][i] + n*tk
    # Solve iteratively kepler equation for eccentricity anomaly
    e = res['Eccentricity'][i]
    Ek = Mk
    counter = 0
    while (counter<100):
        Ek_prev = Ek
        Ek = Ek_prev + (Mk - Ek_prev + e*np.sin(Ek_prev))/(1 - e*np.cos(Ek_prev))
        counter = counter+1
    # Compute true anomaly
    vk = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(Ek/2))
    # Compute argument of latitude 
    fi_k = vk + res['omega'][i]
    # Corrected argument of latitude
    uk = fi_k + res['Cuc'][i]*np.cos(2*fi_k) + res['Cus'][i]*np.sin(2*fi_k)
    # Compute the radial distance rk, considering corrections crc and crs
    rk = A*(1-e*np.cos(Ek)) + res['Crc'][i]*np.cos(2*fi_k) + res['Crs'][i]*np.sin(2*fi_k)
    # Compute the inclination ik of the orbital plane
    ik = res['Io'][i] + res['Cic'][i]*np.cos(2*fi_k) + res['Cis'][i]*np.sin(2*fi_k) + res['IDOT'][i]*tk
    
    # POSITIONS IN OBITAL
    xk_o = rk*np.cos(uk)
    yk_o = rk*np.sin(uk)

    # Compute the longitude of the ascending node λk (with respect to Greenwich)
    omega_k = res['Omega0'][i] + (res['OmegaDot'][i] - omega_dot_E)*tk - omega_dot_E*res['Toe'][i]
    
    # EARTH FIXED POSITIONS
    xk = xk_o*np.cos(omega_k) - yk_o*np.cos(ik)*np.sin(omega_k)
    yk = xk_o*np.sin(omega_k) + yk_o*np.cos(ik)*np.cos(omega_k)
    zk = yk_o*np.sin(ik)
    
    """SATELLITES VELOCITY COMPUTATION"""
    Ek_dot = n/(1-e*np.cos(Ek))
    vk_dot = Ek_dot*np.sqrt(1-e*e)/(1-e*np.cos(Ek))
    ik_dot = res['IDOT'][i] + 2*vk_dot*(res['Cis'][i]*np.cos(2*fi_k) - res['Cic'][i]*np.sin(2*fi_k))
    uk_dot = vk_dot + 2*vk_dot*(res['Cus'][i]*np.cos(2*fi_k) - res['Cuc'][i]*np.sin(2*fi_k))
    rk_dot = e*A*Ek_dot*np.sin(Ek) + 2*vk_dot*(res['Crs'][i]*np.cos(2*fi_k) - res['Crc'][i]*np.sin(2*fi_k))
    omega_k_dot = res['OmegaDot'][i] - omega_dot_E

    # IN PLANE VELOCITIES
    xk_o_dot = rk_dot*np.cos(uk) - rk*uk_dot*np.sin(uk)
    yk_o_dot = rk_dot*np.sin(uk) + rk*uk_dot*np.cos(uk)

    # EARTH FIXED VELOCITIES
    xk_dot = -xk_o*omega_k_dot*np.sin(omega_k) + xk_o_dot*np.cos(omega_k) - yk_o_dot*np.sin(omega_k)*np.cos(ik) - yk_o*(omega_k_dot*np.cos(omega_k)*np.cos(ik) - ik_dot*np.sin(omega_k)*np.sin(ik))
    yk_dot = xk_o*omega_k_dot*np.cos(omega_k) + xk_o_dot*np.sin(omega_k) + yk_o_dot*np.cos(omega_k)*np.cos(ik) - yk_o*(omega_k_dot*np.sin(omega_k)*np.cos(ik) + ik_dot*np.cos(omega_k)*np.sin(ik))
    zk_dot = yk_o_dot*np.sin(ik) + yk_o*ik_dot*np.cos(ik)

    # Inserisco i valori nella tabella dei satelliti
    val = [(res['time'][i], res['sv'][i], xk, yk, zk, xk_dot, yk_dot, zk_dot, dts )]
    new_row = pd.DataFrame(val, columns = ['time', 'sv', 'xs', 'ys', 'zs', 'xs_dot', 'ys_dot', 'zs_dot', 'ts'])
        
    return new_row

def getGPSorbits(nav_path, obs_path, time_range):
    
    nav = gr.load(nav_path, use=['G']).to_dataframe().dropna().sort_values(by=['time'], ascending = True).reset_index()
    nav = fixTime(nav)
    
    obs = gr.load(obs_path, use = 'G', meas = ['C1C'])
    time_obs = obs.time.to_dataframe()['time'].tolist()
    
    satellites = pd.DataFrame(columns = ['time', 'sv', 'xs', 'ys', 'zs', 'xs_dot', 'ys_dot', 'zs_dot', 'ts', 'C1C'])
    
    for w in time_range:
        obs_tk = (((obs.sel(time = w)).to_dataframe()).dropna()).reset_index()
        available_sat = (obs_tk['sv']).to_list()
        if len(available_sat) > 0:
            for sv_i in available_sat:
                # Ogni satellite in vista nell'istante di tempo considerato viene associato allo slot
                # di parametri a cui appartiene attraverso la funzione getPar
                sat_par = getPar(sv_i, w, nav)
                # Se i valori sono disponibli
                if type(sat_par) == pd.core.frame.DataFrame:
                    sat = satPosVel(sat_par).reset_index()
                    sat = sat.merge(obs_tk, on = ['sv', 'time'], how = 'left')
                    satellites = satellites.append(sat)
    satellites = satellites.reset_index().drop(columns=['index'])
    satellites['constellation'] = 'G'
    
    return satellites

def getGALILEOorbits(nav_path, obs_path, time_range):
    nav = gr.load(nav_path, use=['E']).to_dataframe().dropna().sort_values(by=['time'], ascending = True).reset_index()
    nav = fixTime(nav)
    
    tn_in_s = []
    for i in range(len(nav)):
        t = nav['time'][i]
        t_s = t.hour*3600 + t.minute*60 + t.second
        tn_in_s.append(t_s)
        
    nav['time_in_s'] = tn_in_s
    nav = nav[nav['time']>=time_range[0]].reset_index().drop(columns=['index'])
    
    obs = gr.load(obs_path, use = 'E', meas = ['C1C'])
    time_obs = obs.time.to_dataframe()['time'].tolist()
    
    satellites = pd.DataFrame(columns = ['time', 'sv', 'xs', 'ys', 'zs', 'xs_dot', 'ys_dot', 'zs_dot', 'ts', 'C1C'])
    
    for w in time_range:
        w_in_s = w.hour*3600 + w.minute*60 + w.second
        obs_tk = (((obs.sel(time = w)).to_dataframe()).dropna()).reset_index()
        available_sat = (obs_tk['sv']).to_list()
        if len(available_sat) > 0:
            for sv_i in available_sat:
                nav_i = nav[nav['sv'] == sv_i].reset_index().drop(columns=['index'])
                nav_i['delta_T'] = nav_i['time_in_s'] - w_in_s
                nav_i = nav_i[nav_i['delta_T'] >= 0]
                if len(nav_i) >= 0 and nav_i['delta_T'].min()<3600:
                    PAR = nav_i[nav_i['delta_T'] == nav_i['delta_T'].min()].reset_index().drop(columns=['index', 'time'])
                    PAR['time'] = w
                    sat = satPosVel(PAR)
                    sat = sat.merge(obs_tk, on = ['sv', 'time'], how = 'left')
                    satellites = satellites.append(sat)
    
    satellites = satellites.reset_index().drop(columns=['index'])
    satellites['constellation'] = 'E'
    
    return satellites              
        

def checkSatPos(satellites, eph_path):
    
    eph = RINEXreader.readSP3(eph_path)
    eph = eph.rename(columns={'xs':'X', 'ys':'Y', 'zs':'Z', 'ts':'TS'})
        
    '''STATISTICS'''
    
    confronto =  satellites.merge(eph, on = ['sv', 'time']).reset_index()
    confronto['delta_x'] = abs(confronto['xs'] - confronto['X'])
    confronto['delta_y'] = abs(confronto['ys'] - confronto['Y'])
    confronto['delta_z'] = abs(confronto['zs'] - confronto['Z'])
    confronto['delta_ts'] = abs(confronto['ts'] - confronto['TS'])
    confronto = confronto[['time', 'sv', 'delta_x', 'delta_y', 'delta_z', 'delta_ts']]
    
    mean_x = confronto['delta_x'].mean()
    mean_y = confronto['delta_y'].mean()
    mean_z = confronto['delta_z'].mean()
    mean_ts = confronto['delta_ts'].mean()
    x_max = confronto['delta_x'].max()
    y_max = confronto['delta_y'].max()
    z_max =confronto['delta_z'].max()
    ts_max = confronto['delta_ts'].max()
    
    print('abs(computed_value - expected_value for ephemerides')
    print('Mean:')
    print('X: ', mean_x)
    print('Y: ', mean_y)
    print('Z: ', mean_z)
    print('ts: ', mean_ts)
    
    print('Max:')
    print('X: ', x_max)
    print('Y: ', y_max)
    print('Z: ', z_max)
    print('ts: ', ts_max)
    
    return confronto

def checkSatVel(sat, df_gps):
    velocita = pd.DataFrame(columns=['sv', 'time', 'diff_x', 'diff_y', 'diff_z'])
    for i in range(len(sat)):
        row = getPar(sat['sv'][i], sat['time'][i], df_gps)
        row = row.reset_index()
        row = row.drop(columns='index')
        row['time'] = sat['time'][i] + datetime.timedelta(seconds=3)
        check = satPosVel(row.reset_index())
        x_v_t = (check['xs'][0] - sat['xs'][i])/(3)
        y_v_t = (check['ys'][0] - sat['ys'][i])/(3)
        z_v_t = (check['zs'][0] - sat['zs'][i])/(3)
        diff_x = abs(x_v_t - sat['xs_dot'][i])
        diff_y = abs(y_v_t - sat['ys_dot'][i])
        diff_z = abs(z_v_t - sat['zs_dot'][i])
        new_row = pd.DataFrame([[sat['sv'][i], sat['time'][i], diff_x, diff_y, diff_z]], columns=['sv', 'time', 'diff_x', 'diff_y', 'diff_z'])
        velocita = velocita.append(new_row)
    velocita = velocita.reset_index()
    print('abs(computed_value - expected_value for velocities')
    print('Mean:')
    print('X: ', velocita['diff_x'].mean())
    print('Y: ', velocita['diff_y'].mean())
    print('Z: ', velocita['diff_z'].mean())
    
    print('Max:')
    print('X: ', velocita['diff_x'].max())
    print('Y: ', velocita['diff_y'].max())
    print('Z: ', velocita['diff_z'].max())
    return velocita