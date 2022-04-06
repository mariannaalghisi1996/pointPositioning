# -*- coding: utf-8 -*-
"""
@author: Marianna
"""
import pandas as pd
import numpy as np
import math
import codepos.transformations as trf
import codepos.RINEXreader as RINEXreader
import georinex as gr

def ionoCorrectionGPS(phi_u, lambda_u, A, E, GPStime, iono_params):
    alpha = iono_params[0]
    beta = iono_params[1]
    # elevation from 0 to 90 degrees
    E = abs(E)    

    # conversion to semicircles
    phi_u = phi_u / 180
    lambda_u = lambda_u / 180
    A = A / 180
    E = E / 180
    # Speed of light
    c = 2.99792458 * 10**8

    # Psi is the Earth's central angle between the user position and the earth projection of ionospheric intersection point
    psi = 0.0137 / (E + 0.11) - 0.022       
    phi_i = phi_u + psi * math.cos(A * math.pi)

    if abs(phi_i) <= 0.416:
        phi_i = phi_i
    elif phi_i > 0.416:
        phi_i = 0.416
    elif phi_i < -0.416:
        phi_i = -0.416

    # geodetic longitude of the earth projection of the ionospheric intersection point
    lambda_i = lambda_u + psi * math.sin(A * math.pi) / math.cos(phi_i * math.pi)
    # geomagnetic latitude of the earth projection of the ionospheric intersection point
    phi_m = phi_i + 0.064 * math.cos((lambda_i - 1.617) * math.pi)
    # The local time in seconds
    t = 4.32 * 10**4 * lambda_i + GPStime  
    if t >= 86400:
        t = t - 86400
    elif t < 0:
        t = t + 86400

    # Obliquity factor
    F = 1 + 16 * math.pow((0.53 - E), 3)        
    PER = 0
    for n in range(3):
        PER = PER + beta[n] * math.pow(phi_m, n)        

    if PER < 72000:
        PER = 72000    

    x = 2 * math.pi * (t - 50400) / PER
    
    AMP = 0
    for n in range(3):
        AMP  = AMP + alpha[n] * math.pow(phi_m, n)      # the coefficients of a cubic equation representing the amplitude of the vertical delay (4 coefficients - 8 bits each)

    if AMP < 0:
        AMP = 0
    
    if abs(x) >= 1.57:
        T_iono = c *  F * 5 * 10**-9
    else:
        T_iono = c * F * ((5 * 10**-9) + AMP * (1 - (x**2/2) + (x**4 / 24)))

    return T_iono

#def ionoCorrectionGALILEO()

def saastamoinenModel (h, eta):
    if (h > -500) and (h < 5000):
        eta #satellite elevation in radians
        Po = 1013.25 #mBar
        To = 291.15 #degree Kelvin
        Ho = 50/100
        ho = 0
        height = h - ho # h is the ellipsoidal height of the receiver
        Pr = Po * (1-0.0000226 * height)**5.225 #pressure
        Tr = To - 0.0065 * height #temperature
        Hr = Ho * math.exp(-0.0006396 * height)
        er = Hr * math.exp(-37.2465 + (0.213166 * Tr) - 0.000256908 * (Tr)**2) #humidity
        tropoDelay = 0.002277 / math.sin(eta) * (Pr + er * (1255/Tr + 0.05) - (math.tan(eta)) ** -2)
        return tropoDelay
    else:
        tropoDelay = 0
        return tropoDelay

def getTimeCorrections(time_range):
    GPUT = [9.3132257462*10**(-10), 5.329070518*10**(-15), 589824, 2156]
    GLUT = [1.5832483768*10**(-8), 0, 0, 0]
    GPGA = [3.2014213502*10**(-9), -2.220446049*10**(-15), 345600, 2156]
    timeOff = pd.DataFrame(columns = ['time', 'GPST', 'GST', 'GGTO'])
    for w in time_range:
        start_day = w.day - w.weekday() - 1
        t = (w.day - start_day)*24*3600 + w.hour*3600 + w.minute*60 + w.second
        GPST_w = GPUT[0] + GPUT[1]*(t - GPUT[3])
        GST_w = GLUT[0] + GLUT[1]*(t - GLUT[3])
        GGTO_w = GPGA[0] + GPGA[1]*(t - GPGA[3])
        new_row = pd.DataFrame([[w, GPST_w, GST_w, GGTO_w]], columns = ['time', 'GPST', 'GST', 'GGTO'])
        timeOff = timeOff.append(new_row)
    timeOff = timeOff.reset_index().drop(columns=['index'])
    return timeOff
    
def pointPositioning(satellites, nav_path, obs_path, cutoff):
    R_0 = gr.rinexheader(obs_path).get('position')
    iter_max = 20
    omega_dot_E = 7.2921151467*0.00001 #rad/sec
    c = 2.99792458 * 10**8
    
    all_times = []
    for p in satellites['time']:
        if (p not in all_times and p.second == 0):
            all_times.append(p)
    
    timeOff = getTimeCorrections(all_times)
    results_cart = pd.DataFrame(columns = ['datetime', 'xr', 'yr', 'zr', 'dtr_GPS', 'dtr_GAL', 'in_view_sat'])
    
    for w in all_times:
        mask = satellites['time'] == w
        obs_tk = satellites[mask].reset_index().drop(columns=['index'])
        
        start_day = w.day - w.weekday() - 1
        t = (w.day - start_day)*24*3600 + w.hour*3600 + w.minute*60 + w.second
        xr = R_0[0]
        yr = R_0[1]
        zr = R_0[2]
        
        dtr_GPS = 0
        dtr_GAL = 0
        
        for itr in range(iter_max):
            Xr_v = np.array([[xr], [yr], [zr]])
            R0_g = trf.cartToGeod(xr, yr, zr)
            lat0 = trf.degToRad(R0_g[0])
            lon0 = trf.degToRad(R0_g[1])
            h0 = R0_g[2]
            R = np.array([ [-np.sin(lon0), np.cos(lon0), 0],
                [-np.sin(lat0)*np.cos(lon0), -np.sin(lat0)*np.sin(lon0), np.cos(lat0)],
                [np.cos(lat0)*np.cos(lon0), np.cos(lat0)*np.sin(lon0), np.sin(lat0)]
                ])
            azimuth = []
            elevation = []
            for i in range(len(obs_tk)):
                xs = obs_tk['xs'][i]
                ys = obs_tk['ys'][i]
                zs = obs_tk['zs'][i]
                Xs = np.array([[xs],[ys],[zs]])
                delta_X = Xs - Xr_v
                ENU_s = np.dot(R, delta_X)
                E = ENU_s[0][0]
                N = ENU_s[1][0]
                U = ENU_s[2][0]
                hor_dis = np.sqrt(E*E + N*N)
                if hor_dis < 0.1:
                    Az = 0
                    El = 90
                else:
                    Az = trf.radToDeg(np.arctan(E/N))
                    El = trf.radToDeg(np.arctan2(U,hor_dis))
                azimuth.append(Az)
                elevation.append(El)
            obs_tk['Az_S'] = azimuth
            obs_tk['El_S'] = elevation
            
            obs_tk = obs_tk[obs_tk['El_S']>=cutoff].reset_index().drop(columns=['index'])
            
            if len(obs_tk) > 4:
                # Realizzazione matrice A e vettore parametri noti b
                A = []
                iono = []
                tropo = []
                b = []
                sin_el = []
                all_sv = []
                
                for z in range(len(obs_tk)):
                    const_type = obs_tk['constellation'][z]
                    xs = obs_tk['xs'][z]
                    ys = obs_tk['ys'][z]
                    zs = obs_tk['zs'][z]
                    xs_dot = obs_tk['xs_dot'][z]
                    ys_dot = obs_tk['ys_dot'][z]
                    zs_dot = obs_tk['zs_dot'][z]
                    ts = obs_tk['ts'][z]
                    Xs_v = np.array([[xs], [ys], [zs]])
                    Xs_dot_v = np.array([[xs_dot], [ys_dot], [zs_dot]])
                    ro_double = (xs - xr)*(xs - xr) + (ys - yr)*(ys - yr) + (zs - zr)*(zs - zr)
                    ro_approx = np.sqrt(ro_double)
                    e_rs = (Xr_v - Xs_v)/ro_approx
                    prod_vett = Xs_dot_v + np.array([[-omega_dot_E*ys], [omega_dot_E*xs], [0]])
                    alpha_1 = -(e_rs[0][0]*prod_vett[0][0] +  e_rs[1][0]*prod_vett[1][0] + e_rs[2][0]*prod_vett[2][0])
                    alpha_2 = -(e_rs[0][0]*Xs_dot_v[0][0] + e_rs[1][0]*Xs_dot_v[1][0] + e_rs[2][0]*Xs_dot_v[2][0])
                    k_1 = 1 - alpha_1/c + alpha_1*alpha_1/(c*c)
                    k_2 = -alpha_2 + alpha_1*alpha_2/c
                    # effetto relativistico
                    eff_R = 2*(xs*xs_dot + ys*ys_dot + zs*zs_dot)/c
                    
                    # Iono correction
                    iono_params = RINEXreader.getIonoParams(nav_path, 'G')
                    if const_type == 'G':
                        iono_j = ionoCorrectionGPS(trf.radToDeg(lat0), trf.radToDeg(lon0), obs_tk['Az_S'][z], obs_tk['El_S'][z], t, iono_params)
                    else:
                        #iono_j = ionoCorrectionGPS(trf.radToDeg(lat0), trf.radToDeg(lon0), obs_tk['Az_S'][z], obs_tk['El_S'][z], t, iono_params)

                        if 'iono_delay' in obs_tk.columns:
                            iono_j = obs_tk['iono_delay'][z]
                        else:
                            iono_j = 0
                    el_j = trf.degToRad(obs_tk['El_S'][z])
                    tropo_j = saastamoinenModel(h0, el_j)
                    # calcolo b per il j-esimo satellite
                    b_j = k_1*ro_approx - c*ts + eff_R + iono_j + tropo_j
                    b.append([b_j])
                    if const_type == 'G':
                        A.append([k_1*e_rs[0][0], k_1*e_rs[1][0], k_1*e_rs[2][0], (k_2 + c)/c , 0])
                    else:
                        A.append([k_1*e_rs[0][0], k_1*e_rs[1][0], k_1*e_rs[2][0], 0, (k_2 + c)/c ])
                    iono.append(iono_j)
                    tropo.append(tropo_j)
                    all_sv.append(obs_tk['sv'][z])
                    sin_el.append([1/np.sin(el_j)])
                    
                    
                # Termini noti
                P1 = (np.array([obs_tk['C1C'].to_list()])).transpose()
                b = np.array(b)
                A = np.array(A)
                sin_el = np.array(sin_el)
                dP1_oss = P1 - b
                Q = np.identity(len(obs_tk))*sin_el
                Q_inv = np.linalg.inv(Q)
                A_t = A.transpose()
                N = np.dot(np.dot(A_t, Q_inv), A)
                N_inv = np.linalg.inv(N)
                
                # Stima delle incognite
                incognite_stima = np.dot(np.dot(np.dot( N_inv, A_t), Q_inv), dP1_oss )
                
                xr = incognite_stima[0][0] + xr
                yr = incognite_stima[1][0] + yr
                zr = incognite_stima[2][0] + zr
                dtr_GPS = incognite_stima[3][0]/c
                dtr_GAL = incognite_stima[4][0]/c
            
            else:
                xr = np.nan
                yr = np.nan
                zr = np.nan
                dtr_GPS = np.nan
                dtr_GAL = np.nan
            
        new_row = pd.DataFrame([[w, xr, yr, zr, dtr_GPS, dtr_GAL, len(obs_tk)]], columns=['datetime', 'xr', 'yr', 'zr', 'dtr_GPS', 'dtr_GAL', 'in_view_sat'])
        #print(w, 'ok')
        results_cart = results_cart.append(new_row)

    results_cart = results_cart.reset_index()
    return results_cart
