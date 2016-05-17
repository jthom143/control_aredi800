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

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_ohc


## Load Data ###
'''
    Files loaded in this script were created in 'plot_aredi_climatologies.py'
    '''

# Load Temp Data
PATHS={'aredi_400':'~/control_aredi800/data/newCO2_control_400_ag/',
       'aredi_800':'~/control_aredi800/data/newCO2_control_800/',
       'aredi_2400':'~/control_aredi800/data/newCO2_control_2400_ag/',
       'low_gm': '~/control_aredi800/data/gm_min_600/'
        }

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
    if name == 'aredi_2400':
        temp[name] = temp[name][:500]
        rhodzt[name] = rhodzt[name][:500]
    if name == 'low_gm':
        temp[name] = temp[name][:500]
        rhodzt[name] = rhodzt[name][:500]

## Calculate OHC ###
names = {'aredi_400', 'aredi_800', 'aredi_2400', 'low_gm'}

heat = {}
heat_sumz = {}
heat_global = {}

for name in names:
    heat[name], heat_sumz[name], heat_global[name] = calc_ohc(temp[name], rhodzt[name], area[name])

# Plot the variance for the aredi=800 simulation
variance = heat_sumz['aredi_800'].collapsed('time', iris.analysis.VARIANCE).data
lats = heat_sumz['aredi_800'].coord('latitude').points
lons = heat_sumz['aredi_800'].coord('longitude').points


fig1 = plt.figure(figsize = (8, 4))
clevs = np.arange(0, 4.25, 0.25)
#ax = fig1.add_axes([0.1,0.1,0.8,0.9], projection=ccrs.PlateCarree())
plt.contourf(lons, lats, variance/1e19, clevs)
plt.colorbar()
plt.title('A$_{redi}$ = 800 Heat Content Variance')
#plt.savefig('notes/figures/aredi800_heat_variance.png', bbox_inches='tight')

# Plot the variance for all the simulations
fig1 = plt.figure(figsize = (12, 6))
clevs = np.arange(0, 4.25, 0.25)

ax1 = plt.subplot(2,2,1)
variance = heat_sumz['aredi_400'].collapsed('time', iris.analysis.VARIANCE).data
lats = heat_sumz['aredi_400'].coord('latitude').points
lons = heat_sumz['aredi_400'].coord('longitude').points
plt.contourf(lons, lats, variance/1e19, clevs, extend = 'max')
plt.title('A$_{redi}$ = 400')

ax2 = plt.subplot(2,2,2)
variance = heat_sumz['aredi_800'].collapsed('time', iris.analysis.VARIANCE).data
lats = heat_sumz['aredi_800'].coord('latitude').points
lons = heat_sumz['aredi_800'].coord('longitude').points
plt.contourf(lons, lats, variance/1e19, clevs)
plt.title('A$_{redi}$ = 800')

ax3 = plt.subplot(2,2,3)
variance = heat_sumz['aredi_2400'].collapsed('time', iris.analysis.VARIANCE).data
lats = heat_sumz['aredi_2400'].coord('latitude').points
lons = heat_sumz['aredi_2400'].coord('longitude').points
plt.contourf(lons, lats, variance/1e19, clevs)
plt.title('A$_{redi}$ = 2400')

ax4 = plt.subplot(2,2,4)
variance = heat_sumz['low_gm'].collapsed('time', iris.analysis.VARIANCE).data
lats = heat_sumz['low_gm'].coord('latitude').points
lons = heat_sumz['low_gm'].coord('longitude').points
plt.contourf(lons, lats, variance/1e19, clevs)
plt.title('GM$_{min}$ = 600')
plt.savefig('notes/figures/heat_variance.png', bbox_inches='tight')
plt.show()


