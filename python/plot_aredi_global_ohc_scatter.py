"""
    Script to plot ohc scatter plot for different aredi simulations
"""

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib
from netCDF4 import Dataset
from scipy import signal

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_ohc
from calc_convective_periods import convect

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
mld = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
        print name
        temp[name] = iris.load_cube(PATH+'temp.nc')
        rhodzt[name] = iris.load_cube(PATH+'rho_dzt.nc')
        area[name] = iris.load_cube(PATH+'area_t.nc')
	mld[name] = iris.load_cube(PATH+'mld.nc')

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
heat_SH_detrend = {}

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
    heat_SH_detrend[name] = signal.detrend(heat_SH_anomaly[name].data)

## Calculate convective years
year_convect = {}

for name in names:
	year_convect[name], mld_convect, area_convect = convect(mld[name], area[name])

# Figure
fig1 = plt.figure(figsize = (5,10))
ax1 = fig1.add_axes([0.15, 0.75, 0.75, 0.15])
ax1.scatter((global_heat_detrend['aredi_400'])/1e22, heat_SH_detrend['aredi_400']/1e22, marker ='.', color = 'k') 
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11) 

ax2 = fig1.add_axes([0.15, 0.5, 0.75, 0.15])
ax2.scatter((global_heat_detrend['aredi_800'])/1e22, heat_SH_detrend['aredi_800']/1e22, marker ='.', color = 'k')
ax2.set_ylabel('SH Heat Content Anomaly [10$^{22}$ J]')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)

ax3 = fig1.add_axes([0.15, 0.25, 0.75, 0.15])
ax3.scatter((global_heat_detrend['aredi_2400'])/1e22, heat_SH_detrend['aredi_2400']/1e22, marker ='.', color = 'k')
ax3.set_xlabel('Global Heat Content Anomaly [10$^{22}$ J]')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
plt.savefig('notes/aredi_ohc_scatter.png', bbox_inches='tight')



# Figure - highlight convective times
a_400 = ma.masked_where(ma.getmask(year_convect['aredi_400']), global_heat_detrend['aredi_400'])
b_400 = ma.masked_where(ma.getmask(year_convect['aredi_400']), heat_SH_detrend['aredi_400'])

a_800 = ma.masked_where(ma.getmask(year_convect['aredi_800']), global_heat_detrend['aredi_800'])
b_800 = ma.masked_where(ma.getmask(year_convect['aredi_800']), heat_SH_detrend['aredi_800'])

a_2400 = ma.masked_where(ma.getmask(year_convect['aredi_2400']), global_heat_detrend['aredi_2400'])
b_2400 = ma.masked_where(ma.getmask(year_convect['aredi_2400']), heat_SH_detrend['aredi_2400'])

fig1 = plt.figure(figsize = (5,10))
ax1 = fig1.add_axes([0.15, 0.75, 0.75, 0.15])
ax1.scatter((global_heat_detrend['aredi_400'])/1e22, heat_SH_detrend['aredi_400']/1e22, marker ='.', color = 'k')
ax1.scatter(a_400/1e22, b_400/1e22, marker = '.', color = '#368BC1')
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11)

ax2 = fig1.add_axes([0.15, 0.5, 0.75, 0.15])
ax2.scatter((global_heat_detrend['aredi_800'])/1e22, heat_SH_detrend['aredi_800']/1e22, marker ='.', color = 'k')
ax2.scatter(a_800/1e22, b_800/1e22, marker = '.', color = '#368BC1')
ax2.set_ylabel('SH Heat Content Anomaly [10$^{22}$ J]')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)

ax3 = fig1.add_axes([0.15, 0.25, 0.75, 0.15])
ax3.scatter((global_heat_detrend['aredi_2400'])/1e22, heat_SH_detrend['aredi_2400']/1e22, marker ='.', color = 'k')
ax3.scatter(a_2400/1e22, b_2400/1e22, marker = '.', color = '#368BC1')
ax3.set_xlabel('Global Heat Content Anomaly [10$^{22}$ J]')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
plt.savefig('notes/aredi_ohc_scatter_wconv.png', bbox_inches='tight')





plt.show()

