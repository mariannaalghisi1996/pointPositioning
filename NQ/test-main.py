# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 11:52:22 2022

@author: maria
"""
import modelParams as mp
from constants import *
from global_functions import *
import transformations as trf
import TEC as modTEC

import math
import numpy as np
import datetime as datetime
from datetime import datetime as dtt
import pandas as pd

mose_cart = [4642432.50, 1028629.40, 4236854.20]
mose_geod = trf.cartToGeod(4642432.50, 1028629.40, 4236854.20)
lat = mose_geod[0]
lon = mose_geod[1]
h = mose_geod[2]*10**(-3)

time = dtt.strptime('2021-05-06 12:00:00', '%Y-%m-%d %H:%M:%S')

nav_path = 'C:/git/pointPositioning/pointPositioning/mose_nav.rnx'

satellites = pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/sat_galileo.csv')
conv_t = []
for i in range(len(satellites)):
    t_i = satellites['time'][i]
    t_i_conv = dtt.strptime(t_i, '%Y-%m-%d %H:%M:%S')
    conv_t.append(t_i_conv)

satellites = satellites.drop(columns=['time', 'Unnamed: 0'])
satellites['time'] = conv_t

#satellites = satellites[satellites['time'] == time].reset_index().drop(columns=['index'])

P1 = [lat, lon, h]

delays = pd.DataFrame(columns = ['sv', 'time', 'TEC'])
for i in range(len(satellites)):
    time = satellites['time'][i]
    parametri = mp.getModelParams(lat, lon, time, nav_path)
    xs, ys, zs = satellites['xs'][i], satellites['ys'][i], satellites['zs'][i]
    P2 = trf.cartToGeod(xs, ys, zs)
    P2[2] = P2[2]*10**(-3)
    TEC_i = modTEC.getTEC(P1, P2, parametri)
    new_row = pd.DataFrame([[time, satellites['sv'][i], TEC_i]], columns = ['sv', 'time', 'TEC'])
    delays = delays.append(new_row)

delays = delays.reset_index().drop(columns=['index'])

satellites['TEC_new'] = delays['TEC']

satellites.to_csv('C:/git/pointPositioning/pointPositioning/csv/sat_galileo_new.csv')
