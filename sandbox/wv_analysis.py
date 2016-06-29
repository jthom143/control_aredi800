"""
    Script to calculate and plot the wvreenhouse effect
    """

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import matplotlib
from netCDF4 import Dataset
from scipy import signal
import iris.analysis.cartography
import numpy.ma as ma


import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_region import calc_zonal_regions
import colormaps as cmaps

## Load Data ###
'''
    Files loaded in this script were created in 'plot_aredi_climatolowvies.py'
    '''

# Load Temp Data
PATHS={'aredi_800':'~/control_aredi800/data/newCO2_control_800/'}

wv = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    wv[name] = iris.load_cube(PATH+'WVP.nc')

# Define convective and non-convective years based on wvlobal carbon content mins and maxes
p1_min = 136
p1_max = 99

p2_min = 336
p2_max = 204

p5_min = 405
p5_max = 380

wv_p1_contour = wv['aredi_800'][p1_min,:,:] - wv['aredi_800'][p1_max,:,:]
wv_p2_contour = wv['aredi_800'][p2_min,:,:] - wv['aredi_800'][p2_max,:,:]
wv_p5_contour = wv['aredi_800'][p5_min,:,:] - wv['aredi_800'][p5_max,:,:]

wv_mean = wv['aredi_800'].collapsed(['time'], iris.analysis.MEAN)
wv_std = wv['aredi_800'].collapsed(['time'], iris.analysis.STD_DEV)


wv_p1_convection_anomolies = wv['aredi_800'][p1_min,:,:] - wv_mean

lats = wv['aredi_800'].coord('latitude').points
lons = wv['aredi_800'].coord('longitude').points


fig = plt.figure()
clevs = np.arange(0, 51, 1)
ax1 = plt.axes(projection=ccrs.PlateCarree())
im = ax1.contourf(lons, lats, wv_mean.data, clevs, cmap = cmaps.viridis, extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('kg m$^{-2}$')
ax1.coastlines()
plt.title('Water Vapor Climatology', fontsize = 11)

fig = plt.figure()
clevs = np.arange(0, 5.2, 0.1)
ax1 = plt.axes(projection=ccrs.PlateCarree())
im = ax1.contourf(lons, lats, wv_std.data, clevs, cmap = cmaps.viridis, extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('kg m$^{-2}$')
ax1.coastlines()
plt.title('Water Vapor Standard Deviation', fontsize = 11)



fig = plt.figure()
clevs = np.arange(-5, 5.25, 0.25)
ax1 = plt.axes(projection=ccrs.PlateCarree())
im = ax1.contourf(lons, lats, wv_p1_convection_anomolies.data, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('kg m$^{-2}$')
ax1.coastlines()
plt.title('(a) Period 1 Convection', fontsize = 11)

## Mask when not significant (significant when > 2 std)
mask = np.zeros((60,96))
mask2 = np.zeros((60,96))
for i in range(0,len(lats)):
    for j in range(0, len(lons)):
        mask[i,j] = ma.masked_less(wv_p1_convection_anomolies[i,j].data, wv_std[i,j].data)
        mask2[i,j] = ma.masked_greater(wv_p1_convection_anomolies[i,j].data, wv_std[i,j].data*-1)

fig = plt.figure()
clevs = np.arange(-5, 5.25, 0.25)
ax1 = plt.axes(projection=ccrs.PlateCarree())
im = ax1.contourf(lons, lats, wv_p1_convection_anomolies.data, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('kg m$^{-2}$')
ax1.coastlines()
plt.title('(a) Period 1 Convection', fontsize = 11)
im = ax1.contourf(lons, lats, mask, clevs, hatches=['...'], cmap=plt.get_cmap('gray'),
                  extend='both', alpha=0)
im = ax1.contourf(lons, lats, mask2, clevs, hatches=['...'], cmap=plt.get_cmap('gray'),
                  extend='both', alpha=0)
ax1.coastlines()
plt.title('(a) Period 1 Convection', fontsize = 11)


'''
levels = [ 12, wv_p1_contour.data.max() ]
cs = plt.contourf(lons, lats, wv['aredi_800'][p1_min,:,:].data-wv_mean.data, levels = levels, hatches=[None, '.'],
                  cmap=plt.get_cmap('gray'),
                  extend='both', alpha=0)
levels = [ wv_p1_contour.data.max(), -2 ]
cs = plt.contourf(lons, lats, wv_p1_contour.data, levels = levels, hatches=[None, '.'],
                  cmap=plt.get_cmap('gray'),
                  extend='both', alpha=0)
'''
plt.show()
