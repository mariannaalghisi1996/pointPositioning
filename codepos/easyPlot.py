# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 11:19:11 2022

@author: maria
"""

import matplotlib.pyplot as plt
import datetime as datetime
from datetime import datetime as dtt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter
import matplotlib.patches as mpatches

def getPlot(res, x, y, col):
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(res[x], res[y], '-',
            color=col)
    ax.set(xlabel=x, ylabel=y, title = 'Time '+y)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
    plt.show()
    

def getDoublePlot(df1, df2, x, y):
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(df1[x], df1[y], '-', color='blue')
    ax.plot(df2[x], df2[y], '-', color='orange')
    ax.set(xlabel=x, ylabel=y, title = 'Time'+y)
    o_patch = mpatches.Patch(color='blue', label='Klob')
    p_patch = mpatches.Patch(color='orange', label='Klob+Neq')
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
    plt.legend(handles=[p_patch, o_patch])
    plt.show()
    
def getTriplePlot(df1, df2, df3, x, y):
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(df1[x], df1[y], '-', color='orange')
    ax.plot(df2[x], df2[y], '-', color='purple')
    ax.plot(df3[x], df3[y], '-', color='magenta')
    ax.set(xlabel=x, ylabel=y, title = 'Time'+y)
    o_patch = mpatches.Patch(color='orange', label='NeQuickG')
    p_patch = mpatches.Patch(color='purple', label='Klobuchar')
    l_patch = mpatches.Patch(color='magenta', label='iono error = 0')
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
    plt.legend(handles=[p_patch, o_patch, l_patch])
    plt.show()
    
def getTimePlot(res, x, y, z):
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(res[x], res[y], '-', color='pink')
    ax.plot(res[x], res[z], '-', color='purple')
    ax.set(xlabel='time', ylabel='dtr', title = 'Time - dtr')
    o_patch = mpatches.Patch(color='pink', label=y)
    p_patch = mpatches.Patch(color='purple', label=z)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(DateFormatter("%H-%M"))
    plt.legend(handles=[p_patch, o_patch])
    plt.show()