"""
    Script to plot occ variance
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
from calc_OHC_OCC import calc_occ


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

dic = {}
rhodzt = {}
area = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    dic[name] = iris.load_cube(PATH+'dic.nc')
    rhodzt[name] = iris.load_cube(PATH+'rho_dzt.nc')
    area[name] = iris.load_cube(PATH+'area_t.nc')
    
    if name == 'aredi_800':
        dic[name] = dic[name][-500:]
    if name == 'aredi_2400':
        dic[name] = dic[name][:500]
        rhodzt[name] = rhodzt[name][:500]
    if name == 'low_gm':
        dic[name] = dic[name][:500]
        rhodzt[name] = rhodzt[name][:500]

## Calculate OCC ###
names = {'aredi_400', 'aredi_800', 'aredi_2400', 'low_gm'}

carbon = {}
carbon_sumz = {}
carbon_global = {}

for name in names:
    carbon[name], carbon_sumz[name], carbon_global[name] = calc_occ(dic[name], rhodzt[name], area[name])

# Plot the variance for the aredi=800 simulation
variance = carbon_sumz['aredi_800'].collapsed('time', iris.analysis.VARIANCE).data
lats = carbon_sumz['aredi_800'].coord('latitude').points
lons = carbon_sumz['aredi_800'].coord('longitude').points


fig1 = plt.figure(figsize = (8, 4))
clevs = np.arange(0, 105, 1)
#ax = fig1.add_axes([0.1,0.1,0.8,0.9], projection=ccrs.PlateCarree())
plt.contourf(lons, lats, variance/1e3, clevs)
plt.colorbar()
plt.title('A$_{redi}$ = 800 Carbon Content Variance')
#plt.savefig('notes/figures/aredi800_carbon_variance.png', bbox_inches='tight')

# Plot the variance for all the simulations
fig1 = plt.figure(figsize = (12, 6))
clevs = np.arange(0, 105, 1)

ax1 = plt.subplot(2,2,1)
variance = carbon_sumz['aredi_400'].collapsed('time', iris.analysis.VARIANCE).data
lats = carbon_sumz['aredi_400'].coord('latitude').points
lons = carbon_sumz['aredi_400'].coord('longitude').points
plt.contourf(lons, lats, variance/1e3, clevs)
plt.title('A$_{redi}$ = 400')

ax2 = plt.subplot(2,2,2)
variance = carbon_sumz['aredi_800'].collapsed('time', iris.analysis.VARIANCE).data
lats = carbon_sumz['aredi_800'].coord('latitude').points
lons = carbon_sumz['aredi_800'].coord('longitude').points
plt.contourf(lons, lats, variance/1e3, clevs)
plt.title('A$_{redi}$ = 800')

ax3 = plt.subplot(2,2,3)
variance = carbon_sumz['aredi_2400'].collapsed('time', iris.analysis.VARIANCE).data
lats = carbon_sumz['aredi_2400'].coord('latitude').points
lons = carbon_sumz['aredi_2400'].coord('longitude').points
plt.contourf(lons, lats, variance/1e3, clevs)
plt.title('A$_{redi}$ = 2400')

ax4 = plt.subplot(2,2,4)
variance = carbon_sumz['low_gm'].collapsed('time', iris.analysis.VARIANCE).data
lats = carbon_sumz['low_gm'].coord('latitude').points
lons = carbon_sumz['low_gm'].coord('longitude').points
plt.contourf(lons, lats, variance/1e3, clevs)
plt.title('GM$_{min}$ = 600')
plt.savefig('notes/figures/carbon_variance.png', bbox_inches='tight')
plt.show()


