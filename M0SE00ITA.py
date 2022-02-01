# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 09:35:18 2022

@author: Marianna

Code for the computation of the geodetic coordinates of Mose at a given epoch.
Mose belongs to EPN, through EPN website it is possible to obtain Mose's position 
and velocities relative tu ETRF14. The proposed values refers to a specific epoch
t_0 = 2010.0 and, in particular, the proposed coordinates at t_0 are different 
with respect the coordinates at epoch t global geodynamics and local deformations. 
"""

import datetime as datetime
from datetime import datetime as dtt

'''ETRF2014 reference values'''
X_0 = 4642432.712
Y_0 = 1028629.121
Z_0 = 4236854.055

v_X = -0.0009
v_Y = -0.0014
v_Z = 0.0003

# Starting epoch
t_0 = 2010.0

#Current epoch
time = dtt.strptime('2021-05-06', '%Y-%m-%d')
t = time.year + ((time.month-1)*30 + time.day)/365

delta_t = t - t_0

X = X_0 + v_X*delta_t
Y = Y_0 + v_Y*delta_t
Z = Z_0 + v_Z*delta_t

'''
Final results for t = 2021-05-06:
    X = 4642432.701789316
    Y = 1028629.1051167124
    Z = 4236854.058403561
'''