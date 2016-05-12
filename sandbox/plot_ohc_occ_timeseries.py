"""
    Script to plot the newCO2 control annual ocean heat content and ocean carbon content
"""

import iris
import matplotlib.pyplot as plt
import numpy as np
import numpy as np
from scipy import signal
import iris.analysis.cartography

# Import my functions
import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_ohc, calc_occ
from calc_region import calc_weddell_sea

# Data Path:
PATH = '~/control_aredi800/data/newCO2_control_800/'

# Load Data:
dic = iris.load_cube(PATH+'dic.nc')
temp = iris.load_cube(PATH+'temp.nc')
rhodzt = iris.load_cube(PATH+'rho_dzt.nc')
area = iris.load_cube(PATH+'area_t.nc')
mld = iris.load_cube(PATH+'mld.nc')

# Ocean Bling has longer timeseries, trim dic down to same length as temp
dic = dic[140:,:,:,:]
if dic[0,:,:,:].coord('time').points != temp[0,:,:,:].coord('time').points:
    raise Exception('dic and temp do not have the same initial time')

# Calculate OHC and OCC
heat, heat_sumz, heat_global = calc_ohc(temp, rhodzt, area)
dic, dic_sumz, dic_global = calc_occ(dic, rhodzt, area)

# Calculate temp in Weddell Sea at depth: 1500-2500m
temp_ws = calc_weddell_sea(temp)
temp_ws = temp_ws[:,18:22,:,:]
temp_ws.coord('latitude').guess_bounds()
temp_ws.coord('longitude').guess_bounds()
temp_ws.coord('tcell pstar').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(temp_ws)
temp_convect = temp_ws.collapsed(['tcell pstar', 'latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)

# Calculate mld in Weddell Sea
mld_ws = calc_weddell_sea(mld)
mld_ws.coord('latitude').guess_bounds()
mld_ws.coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(mld_ws)
mld_ws = mld_ws.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)

# Southern Ocean
SH = iris.Constraint(latitude=lambda y: -90 < y < -50)
dic_sumz_SH = dic_sumz.extract(SH)
heat_sumz_SH = heat_sumz.extract(SH)
area_SH = area.extract(SH)

# Sum over Southern Ocean
dic_sumz_weighted = dic_sumz_SH*area_SH
dic_SH_total = dic_sumz_weighted.collapsed(['longitude', 'latitude'], iris.analysis.SUM)

heat_sumz_weighted = heat_sumz_SH*area_SH
heat_SH_total = heat_sumz_weighted.collapsed(['longitude', 'latitude'], iris.analysis.SUM)

# Time mean
dic_SH_mean = dic_SH_total.collapsed('time', iris.analysis.MEAN).data
heat_SH_mean = heat_SH_total.collapsed('time', iris.analysis.MEAN).data

heat_mean = heat_global.collapsed('time', iris.analysis.MEAN).data
dic_mean = dic_global.collapsed('time', iris.analysis.MEAN).data

# Detrend
#Detrend global
dic_anomaly = dic_global.data-dic_mean
dic_anomaly_detrend = signal.detrend(dic_anomaly)

dic_SH_anomaly = dic_SH_total.data-dic_SH_mean
dic_SH_anomaly_detrend = signal.detrend(dic_SH_anomaly)


# Timeseries
fig = plt.figure(figsize = (8, 10))
ax1 = fig.add_axes([0.1,0.7,0.8,0.22])
plt.title('(a) Convection', fontsize = 12)
ax1.plot(mld_ws.data, color = 'k', lw = 1.5)
plt.gca().invert_yaxis()
ax2 = ax1.twinx()
ax2.plot(temp_convect.data, 'g', lw = 1.5)
ax1.set_ylabel('Depth (m)')
ax2.set_ylabel('Degrees C')
plt.xticks([0,100, 200, 300, 400, 500], [' ',' ', ' ', ' ', ' ', ' '])


ax1 = fig.add_axes([0.1,0.42,0.8,0.22])
ax1.plot((heat_SH_total.data-heat_SH_mean)/1e22, color = 'b', lw = 1.5, label = 'Southern Ocean')
ax1.plot((heat_global.data-heat_mean)/1e22, color = 'k', lw = 1.5, label = 'Global')
plt.title('(b) Heat Content', fontsize = 12)
plt.axhline(0, color = 'k', ls = '--')
plt.ylabel('$10^{22}$ Joules')
plt.ylim([-6, 6])
plt.xticks([0,100, 200, 300, 400, 500], [' ',' ', ' ', ' ', ' ', ' '])


ax1 = fig.add_axes([0.1,0.12,0.8,0.22])
ax1.plot((dic_SH_anomaly_detrend)/1e15, color = 'b', lw = 1.5, label = 'Southern Ocean')
ax1.plot(dic_anomaly_detrend/1e15, color = 'k', lw = 1.5, label = 'Global')
plt.title('(c) Carbon Content', fontsize = 12)
plt.ylabel('Pg Carbon')
plt.axhline(0, color = 'k', ls = '--')
plt.legend(fontsize = 11)
plt.ylim([-5, 5])
ax1.set_xlabel('Time (years)')




plt.show()





