"""
    Script to plot ohc variance
    """

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import matplotlib.pyplot as plt
import numpy as np
import numpy as np
import cartopy.crs as ccrs
import matplotlib
from netCDF4 import Dataset
from scipy import signal
from iris.analysis.stats import pearsonr

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_ohc


## Load Data ###
'''
    Files loaded in this script were created in 'plot_aredi_climatologies.py'
    '''

# Load Temp Data
PATHS={'aredi_800':'~/control_aredi800/data/newCO2_control_800/'}

temp = {}
rhodzt = {}
area = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    temp[name] = iris.load_cube(PATH+'temp.nc')
    rhodzt[name] = iris.load_cube(PATH+'rho_dzt.nc')
    area[name] = iris.load_cube(PATH+'area_t.nc')
    
    if name == 'aredi_800':
        temp[name] = temp[name][-500:]

# Isolate the top 1000m
temp['aredi_800'] = temp['aredi_800'][:,:16,:,:]
rhodzt['aredi_800'] = rhodzt['aredi_800'][:,:16,:,:]

## Calculate OHC ###
names = {'aredi_800'}

heat = {}
heat_sumz = {}
heat_global = {}

for name in names:
    heat[name], heat_sumz[name], heat_global[name] = calc_ohc(temp[name], rhodzt[name], area[name])

## Plot ohc maps in the top 1000m for convective - non-convective
p1_min = 136
p1_max = 99

p2_min = 337
p2_max = 204

p5_min = 405
p5_max = 380

heat_p1_contour = heat_sumz['aredi_800'][p1_min,:,:] - heat_sumz['aredi_800'][p1_max,:,:]
heat_p2_contour = heat_sumz['aredi_800'][p2_min,:,:] - heat_sumz['aredi_800'][p2_max,:,:]
heat_p5_contour = heat_sumz['aredi_800'][p5_min,:,:] - heat_sumz['aredi_800'][p5_max,:,:]


lats = heat_sumz['aredi_800'].coord('latitude').points
lons = heat_sumz['aredi_800'].coord('longitude').points


fig = plt.figure()
clevs = np.arange(-6, 6.5, 0.5)
ax1 = plt.axes(projection=ccrs.PlateCarree())
im = ax1.contourf(lons, lats, heat_p1_contour.data/1e9, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('10$^9$ Joules')
ax1.coastlines()
plt.title('(a) Period 1', fontsize = 11)
plt.savefig('notes/figures/heat_top1000_contour_1.png', bbox_inches='tight')



fig = plt.figure()
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.contourf(lons, lats, heat_p1_contour.data/1e9, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('10$^9$ Joules')
ax2.coastlines()
plt.title('(b) Period 2', fontsize = 11)
plt.savefig('notes/figures/heat_top1000_contour_2.png', bbox_inches='tight')


fig = plt.figure()
ax3 = plt.axes(projection=ccrs.PlateCarree())
ax3.contourf(lons, lats, heat_p1_contour.data/1e9, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('10$^9$ Joules')
ax3.coastlines()
plt.title('(c) Period 3', fontsize = 11)
plt.savefig('notes/figures/heat_top1000_contour_3.png', bbox_inches='tight')


### Load in SST contour data
sst_p1_contour = iris.load_cube('data/sst_p1_contour.nc')
sst_p2_contour = iris.load_cube('data/sst_p2_contour.nc')
sst_p5_contour = iris.load_cube('data/sst_p5_contour.nc')

sst = iris.load_cube('data/sst_cube.nc')

heat_sst_corr = pearsonr(heat_sumz['aredi_800'], sst, corr_coords={'time'})

fig = plt.figure()
clevs = np.arange(-1, 1.1, 0.1)
ax1 = plt.axes(projection=ccrs.PlateCarree())
im = ax1.contourf(lons, lats, heat_sst_corr.data, clevs,  cmap = 'RdBu_r', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
plt.title('Correlation between SST and OHC in top 900m', fontsize = 11)
ax1.coastlines()


### Load in g contour data
g_p1_contour = iris.load_cube('data/g_p1_contour.nc')
g_p2_contour = iris.load_cube('data/g_p2_contour.nc')
g_p5_contour = iris.load_cube('data/g_p5_contour.nc')

heat_g_corr = pearsonr(heat_p1_contour, g_p1_contour, corr_coords={'longitude'})

fig = plt.figure()
plt.plot(heat_g_corr.data, lats, color = 'k')
plt.xlim([-1, 1])
plt.ylim([-90, 90])
ax3.axvline(0, color = 'k', ls = '--')






plt.show()