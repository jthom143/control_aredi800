"""
    Script to calculate and plot the greenhouse effect
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

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_region import calc_zonal_regions


## Load Data ###
'''
    Files loaded in this script were created in 'plot_aredi_climatologies.py'
'''

# Load Temp Data
PATHS={'aredi_800':'~/control_aredi800/data/newCO2_control_800/'}

olr = {}
lwup_sfc = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    olr[name] = iris.load_cube(PATH+'olr.nc')
    lwup_sfc[name] = iris.load_cube(PATH+'lwup_sfc.nc')


G = lwup_sfc['aredi_800'] - olr['aredi_800']
G_zonal = G.collapsed('longitude', iris.analysis.MEAN)
G_zonal_mean = G_zonal.collapsed('time', iris.analysis.MEAN)

# Define convective and non-convective years based on global carbon content mins and maxes
p1_min = 136
p1_max = 99

p2_min = 336
p2_max = 204

p5_min = 405
p5_max = 380

G_p1 = G_zonal[p1_min,:] - G_zonal[p1_max,:]
G_p2 = G_zonal[p2_min,:] - G_zonal[p2_max,:]
G_p5 = G_zonal[p5_min,:] - G_zonal[p5_max,:]

G_p1_contour = G[p1_min,:,:] - G[p1_max,:,:]
G_p2_contour = G[p2_min,:,:] - G[p2_max,:,:]
G_p5_contour = G[p5_min,:,:] - G[p5_max,:,:]

lats = olr['aredi_800'].coord('latitude').points
lons = olr['aredi_800'].coord('longitude').points

f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
ax1.plot(lats,G_p1.data, color = 'k', lw = 2)
ax1.axhline(0, color = 'k', ls = '--')

ax2.plot(lats, G_p2.data, color = 'k', lw = 2)
ax2.axhline(0, color = 'k', ls = '--')

ax3.plot(lats, G_p5.data, color = 'k', lw = 2)
ax3.axhline(0, color = 'k', ls = '--')
ax3.set_xlim([-90,90])
plt.savefig('notes/figures/greenhouse_zonal.png', bbox_inches='tight')


f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
ax1.plot(lats, G_zonal[p1_min].data - G_zonal_mean.data, color = 'r', lw = 2)
ax1.plot(lats, G_zonal[p1_max].data - G_zonal_mean.data, color = 'b', lw = 2)
ax1.axhline(0, color = 'k', ls = '--')

ax2.plot(lats, G_zonal[p2_min].data - G_zonal_mean.data, color = 'r', lw = 2)
ax2.plot(lats, G_zonal[p2_max].data - G_zonal_mean.data, color = 'b', lw = 2)
ax2.axhline(0, color = 'k', ls = '--')

ax3.plot(lats, G_zonal[p5_min].data - G_zonal_mean.data, color = 'r', lw = 2)
ax3.plot(lats, G_zonal[p5_max].data - G_zonal_mean.data, color = 'b', lw = 2)
ax3.axhline(0, color = 'k', ls = '--')
ax3.set_xlim([-90,90])

f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex='col', sharey='row')
ax1.plot(lats, G_zonal[p1_min].data, color = 'r', lw = 2)
ax1.plot(lats, G_zonal[p1_max].data , color = 'b', lw = 2)

ax2.plot(lats, G_zonal[p2_min].data, color = 'r', lw = 2)
ax2.plot(lats, G_zonal[p2_max].data , color = 'b', lw = 2)

ax3.plot(lats, G_zonal[p5_min].data, color = 'r', lw = 2)
ax3.plot(lats, G_zonal[p5_max].data, color = 'b', lw = 2)
ax3.set_xlim([-90,90])


fig = plt.figure()
clevs = np.arange(-40, 42, 2)
ax1 = plt.axes(projection=ccrs.PlateCarree())
im = ax1.contourf(lons, lats, G_p1_contour.data, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('W m$^{-2}$')
ax1.coastlines()
plt.title('(a) Period 1', fontsize = 11)
plt.savefig('notes/figures/greenhouse_contour_1.png', bbox_inches='tight')



fig = plt.figure()
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.contourf(lons, lats, G_p2_contour.data, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('W m$^{-2}$')
ax2.coastlines()
plt.title('(b) Period 2', fontsize = 11)
plt.savefig('notes/figures/greenhouse_contour_2.png', bbox_inches='tight')


fig = plt.figure()
ax3 = plt.axes(projection=ccrs.PlateCarree())
ax3.contourf(lons, lats, G_p5_contour.data, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('W m$^{-2}$')
ax3.coastlines()
plt.title('(c) Period 3', fontsize = 11)
plt.savefig('notes/figures/greenhouse_contour_3.png', bbox_inches='tight')


##### Calculate normalized greenhouse effect (g)
SB = iris.load_cube('data/stephan_boltzmann.nc')
# Regrid G onto ocean grid
G_new = G.regrid(SB, iris.analysis.Linear())

g = G_new/SB

g_p1_contour = g[p1_min,:,:] - g[p1_max,:,:]
g_p2_contour = g[p2_min,:,:] - g[p2_max,:,:]
g_p5_contour = g[p5_min,:,:] - g[p5_max,:,:]

lats = g.coord('latitude').points
lons = g.coord('longitude').points


fig = plt.figure()
clevs = np.arange(-0.15, 0.16, 0.01)
ax1 = plt.axes(projection=ccrs.PlateCarree())
im = ax1.contourf(lons, lats, g_p1_contour.data, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('W m$^{-2}$')
ax1.coastlines()
plt.title('(a) Period 1', fontsize = 11)
#plt.savefig('notes/figures/greenhouse_contour_1.png', bbox_inches='tight')



fig = plt.figure()
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.contourf(lons, lats, g_p2_contour.data, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('W m$^{-2}$')
ax2.coastlines()
plt.title('(b) Period 2', fontsize = 11)
#plt.savefig('notes/figures/greenhouse_contour_2.png', bbox_inches='tight')


fig = plt.figure()
ax3 = plt.axes(projection=ccrs.PlateCarree())
ax3.contourf(lons, lats, g_p5_contour.data, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('W m$^{-2}$')
ax3.coastlines()
plt.title('(c) Period 3', fontsize = 11)
#plt.savefig('notes/figures/greenhouse_contour_3.png', bbox_inches='tight')

## Save data
iris.save(g_p1_contour, 'data/g_p1_contour.nc')
iris.save(g_p2_contour, 'data/g_p2_contour.nc')
iris.save(g_p5_contour, 'data/g_p5_contour.nc')

iris.save(g, 'data/g.nc')


plt.show()
