"""
    Script to plot the timeseries of carbon and heat in the weddell sea, southern ocean, and global for the aredi 800 simulation
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
from calc_OHC_OCC import calc_occ, calc_ohc
from calc_region import calc_weddell_sea


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
temp = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    dic[name] = iris.load_cube(PATH+'dic.nc')
    rhodzt[name] = iris.load_cube(PATH+'rho_dzt.nc')
    area[name] = iris.load_cube(PATH+'area_t.nc')
    temp[name] = iris.load_cube(PATH+'temp.nc')
    
    if name == 'aredi_800':
        dic[name] = dic[name][-500:]
    if name == 'aredi_2400':
        dic[name] = dic[name][:500]
        rhodzt[name] = rhodzt[name][:500]
        temp[name] = temp[name][:500]
    if name == 'low_gm':
        dic[name] = dic[name][:500]
        rhodzt[name] = rhodzt[name][:500]
        temp[name] = temp[name][:500]


## Calculate OCC ###
names = {'aredi_400', 'aredi_800', 'aredi_2400', 'low_gm'}

carbon = {}
carbon_sumz = {}
carbon_global = {}
global_carbon_mean = {}
global_carbon_anomaly = {}
global_carbon_detrend = {}
carbon_SH = {}
carbon_SH_anomaly = {}
carbon_SH_detrend = {}

carbon_ws = {}
carbon_ws_anomaly= {}

for name in names:
    carbon[name], carbon_sumz[name], carbon_global[name] = calc_occ(dic[name], rhodzt[name], area[name])
    global_carbon_mean[name] = carbon_global[name].collapsed('time', iris.analysis.MEAN)
    global_carbon_anomaly[name] = carbon_global[name].data-global_carbon_mean[name].data
    
    # Weddell Sea:
    carbon_ws[name]= calc_weddell_sea(carbon_sumz[name])
    carbon_ws[name].coord('latitude').guess_bounds()
    carbon_ws[name].coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(carbon_ws[name])
    carbon_ws[name] = carbon_ws[name].collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights = grid_areas)
    carbon_ws_mean = carbon_ws[name].collapsed('time', iris.analysis.MEAN)
    carbon_ws_anomaly[name] = carbon_ws[name].data - carbon_ws_mean.data

    
    # Southern Hemisphere:
    constraint = iris.Constraint(latitude=lambda y: -90 < y < -50)
    carbon_SH_tmp = carbon_sumz[name].extract(constraint)
    area_SH = area[name].extract(constraint)
    b = carbon_SH_tmp*area_SH
    carbon_SH[name] = b.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    mean_SH = carbon_SH[name].collapsed('time', iris.analysis.MEAN)
    carbon_SH_anomaly[name] = carbon_SH[name] - mean_SH
    
    # Detrend with polynomial
    x = np.arange(0, 500, 1)
    
    p1 = np.poly1d(np.polyfit(x, global_carbon_anomaly[name]/1e15, 4))
    global_carbon_detrend[name] = global_carbon_anomaly[name]/1e15 - p1(x)
    
    p2 = np.poly1d(np.polyfit(x, carbon_SH_anomaly[name].data/1e15, 4))
    carbon_SH_detrend[name] = carbon_SH_anomaly[name].data/1e15 - p2(x)

fig1 = plt.figure(figsize = (8, 4))
plt.plot(global_carbon_detrend['aredi_800'], color = 'k', lw = 2)
plt.plot(carbon_SH_detrend['aredi_800'], color = 'k', lw = 2, ls = '--')
plt.plot(carbon_ws_anomaly['aredi_800']/1e15, color = 'k', lw = 2, ls = '-.')
#plt.savefig('notes/figures/aredi800_occ_timeseries.png')

### Calculate OHC
heat = {}
heat_sumz = {}
heat_global = {}
global_heat_mean = {}
global_heat_anomaly = {}
global_heat_detrend = {}
heat_SH = {}
heat_SH_anomaly = {}
heat_SH_detrend = {}

heat_ws = {}
heat_ws_anomaly= {}

for name in names:
    heat[name], heat_sumz[name], heat_global[name] = calc_ohc(temp[name], rhodzt[name], area[name])
    global_heat_mean[name] = heat_global[name].collapsed('time', iris.analysis.MEAN)
    global_heat_anomaly[name] = heat_global[name].data-global_heat_mean[name].data
    
    # Southern Hemisphere:
    constraint = iris.Constraint(latitude=lambda y: -90 < y < -50)
    heat_SH_tmp = heat_sumz[name].extract(constraint)
    area_SH = area[name].extract(constraint)
    b = heat_SH_tmp*area_SH
    heat_SH[name] = b.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    mean_SH = heat_SH[name].collapsed('time', iris.analysis.MEAN)
    heat_SH_anomaly[name] = heat_SH[name] - mean_SH
    
    # Weddell Sea:
    heat_ws[name]= calc_weddell_sea(heat_sumz[name])
    heat_ws[name].coord('latitude').guess_bounds()
    heat_ws[name].coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(heat_ws[name])
    heat_ws[name] = heat_ws[name].collapsed(['latitude', 'longitude'], iris.analysis.SUM, weights = grid_areas)
    heat_ws_mean = heat_ws[name].collapsed('time', iris.analysis.MEAN)
    heat_ws_anomaly[name] = heat_ws[name].data - heat_ws_mean.data
    
    # Detrend with polynomial
    x = np.arange(0, 500, 1)
    
    p1 = np.poly1d(np.polyfit(x, global_heat_anomaly[name]/1e22, 4))
    global_heat_detrend[name] = global_heat_anomaly[name]/1e22 - p1(x)
    
    p2 = np.poly1d(np.polyfit(x, heat_SH_anomaly[name].data/1e22, 4))
    heat_SH_detrend[name] = heat_SH_anomaly[name].data/1e22 - p2(x)


fig1 = plt.figure(figsize = (8, 4))
plt.plot(global_heat_detrend['aredi_800'], color = 'k', lw = 2)
plt.plot(heat_SH_detrend['aredi_800'], color = 'k', lw = 2, ls = '--')
plt.plot(heat_ws_anomaly['aredi_800']/1e22, color = 'k', lw = 2, ls = '-.')




plt.show()
