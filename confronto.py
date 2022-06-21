# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 11:37:10 2022

@author: maria
"""

import pandas as pd
import codepos.easyPlot as ep
import codepos.transformations as trf

path = 'C:/git/pointPositioning/pointPositioning/testRINEX/'

rLC2 = pd.read_csv(path+'results2_LC.csv')
rLC3 = pd.read_csv(path+'results3_LC.csv')

cLC = pd.DataFrame()
cLC['dE'] = rLC2['E'] - rLC3['E']
cLC['dN'] = rLC2['N'] - rLC3['N']
cLC['dU'] = rLC2['U'] - rLC3['U']

cLC['time'] = rLC2['time']
cLC = trf.fixDateTime(cLC)

ep.getPlot(cLC, 'time', 'dE', 'red')
ep.getPlot(cLC, 'time', 'dN', 'blue')
ep.getPlot(cLC, 'time', 'dU', 'green')




