# -*- coding: utf-8 -*-
"""
@author: Marianna
"""

import georinex as gr
import pandas as pd
import datetime as datetime
from datetime import datetime as dtt
import os

def convertDate(date):
    date = (date)[1:(len(date[22])-2)]
    date = date.split(' ')
    date_list = []
    for i in date:
        if i != '':
            date_list.append(i)
    date_string = date_list[0]+'-'+date_list[1]+'-'+date_list[2]+' '+date_list[3]+':'+date_list[4]+':'+date_list[5].split('.')[0]
    date_time = dtt.strptime(date_string, '%Y-%m-%d %H:%M:%S')
    return date_time

def readSP3(file_path):
    eph = pd.DataFrame(columns=['time', 'sv', 'xs', 'ys', 'zs', 'ts'])
    file = open(file_path, 'r')
    lines = file.readlines()
    for i in range(22, len(lines)):
       line = lines[i]
       if line[0] == '*':
           date = convertDate(line)
       else:
           if (line[0] == 'P'):
               line_list = line.split()
               if line_list[4] != '9999':
                   sv_i = line_list[0][1:]
                   xs_i = float(line_list[1])*1000
                   ys_i = float(line_list[2])*1000
                   zs_i = float(line_list[3])*1000
                   ts_i = float(line_list[4])*0.000001
                   
                   new_eph = pd.DataFrame([[date, sv_i, xs_i, ys_i, zs_i, ts_i]], columns=['time', 'sv', 'xs', 'ys', 'zs', 'ts'])
                   eph = eph.append(new_eph)
    eph = eph.reset_index().drop(columns = ['index'])
    return eph
    
def getIonoParams(file_path):
    file = open(file_path, 'r')
    lines = file.readlines()
    ion_alpha = lines[5].split()[0:4]
    ion_beta = lines[6].split()[0:4]
    for i in range(len(ion_alpha)):
        ion_alpha[i] = float(ion_alpha[i])
        ion_beta[i] = float(ion_beta[i])
    return [ion_alpha, ion_beta]

def getStartPos(file_path):
    file = open(file_path, 'r')
    lines = file.readlines()
    line_pos = lines[9].split()[0:3]
    for i in range(3):
        line_pos[i] = float(line_pos[i])
    return line_pos
          