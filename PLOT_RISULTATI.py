# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 09:28:19 2022

@author: maria
"""

import folium
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

import transformations as trf

sat_gal = pd.read_csv('C:/git/pointPositioning/pointPositioning/csv/sat_galileo_new.csv')

lat, lon, h = [], [], []

for i in range(len(sat_gal)):
    geod = trf.cartToGeod(sat_gal['xs'][i], sat_gal['ys'][i], sat_gal['xs'][i])
    lat.append(geod[0])
    lon.append(geod[1])
    h.append(geod[2])

sat_gal['Latitude'] = lat
sat_gal['Longitude'] = lon
sat_gal['h'] = h

sat_E3 = sat_gal[sat_gal['sv'] == 'E03'].reset_index()

M0SE_cart = [4642432.701789316, 1028629.1051167124, 4236854.058403561]
M0SE_geod = trf.cartToGeod(4642432.701789316, 1028629.1051167124, 4236854.058403561)

geometry = [Point(xy) for xy in zip(sat_E3['Longitude'], sat_E3['Latitude'])]
gdf = gpd.GeoDataFrame(sat_E3, geometry=geometry)

gdf.to_file('C:/git/pointPositioning/pointPositioning/shapefile/E03.shp')


