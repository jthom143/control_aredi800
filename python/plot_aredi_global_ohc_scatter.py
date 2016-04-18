"""
    Script to plot ohc scatter plot for different aredi simulations
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
       'aredi_2400':'~/control_aredi800/data/newCO2_control_2400_ag/'
              #'low_gm': 'data/derived/ocean_temp_low_gm.nc'
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


## Calculate OHC ###
names = {'aredi_400', 'aredi_800', 'aredi_2400'}

heat = {}
heat_sumz = {}
heat_global = {}
global_heat_mean = {}
global_heat_anomaly = {}
global_heat_detrend = {}
heat_SH = {}
heat_SH_anomaly = {}

for name in names:
    heat[name], heat_sumz[name], heat_global[name] = calc_ohc(temp[name], rhodzt[name], area[name])
    global_heat_mean[name] = heat_global[name].collapsed('time', iris.analysis.MEAN)
    global_heat_anomaly[name] = heat_global[name].data-global_heat_mean[name].data
    global_heat_detrend[name] = signal.detrend(global_heat_anomaly[name])

    # Southern Hemisphere:
    constraint = iris.Constraint(latitude=lambda y: -90 < y < -50)
    heat_SH_tmp = heat_sumz[name].extract(constraint)
    area_SH = area[name].extract(constraint)
    b = heat_SH_tmp*area_SH
    heat_SH[name] = b.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    mean_SH = heat_SH[name].collapsed('time', iris.analysis.MEAN)
    heat_SH_anomaly[name] = heat_SH[name] - mean_SH

# Figure
fig1 = plt.figure(figsize = (5,10))
ax1 = fig1.add_axes([0.15, 0.75, 0.75, 0.15])
ax1.scatter((global_heat_detrend['aredi_400'])/1e22, heat_SH_anomaly['aredi_400'].data/1e22, marker ='.', color = 'k') 
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11) 

ax2 = fig1.add_axes([0.15, 0.5, 0.75, 0.15])
ax2.scatter((global_heat_detrend['aredi_800'])/1e22, heat_SH_anomaly['aredi_800'].data/1e22, marker ='.', color = 'k')
ax2.set_ylabel('SH Heat Content Anomaly [10$^{22}$ J]')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)

ax3 = fig1.add_axes([0.15, 0.25, 0.75, 0.15])
ax3.scatter((global_heat_detrend['aredi_2400'])/1e22, heat_SH_anomaly['aredi_2400'].data/1e22, marker ='.', color = 'k')
ax3.set_xlabel('Global Heat Content Anomaly [10$^{22}$ J]')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
plt.savefig('notes/aredi_ohc_scatter.png', bbox_inches='tight')
plt.show()
