"""
    Script to plot the newCO2 control annual ocean heat content and ocean carbon content
"""

import iris
import matplotlib.pyplot as plt
import numpy as np
import numpy as np
from scipy import signal

# Import my functions
import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_OHC, calc_OCC

# Data Path:
PATH = '~/control_aredi800/data/'

# Load Data:
dic = iris.load_cube(PATH+'dic.nc')
temp = iris.load_cube(PATH+'temp.nc')
rhodzt = iris.load_cube(PATH+'rhodzt')
area = iris.load_cube(PATH+'area_t.nc')

# Calculate OHC and OCC
heat, heat_sumz, heat_global = calc_OHC(temp, rhodzt, area)
dic, dic_sumz, dic_global = calc_OCC(dic, rhodzt, area)

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
dic_anomaly_detrend = signal.detrend(dic_anomoly)

dic_SH_anomaly = dic_SH_total.data-dic_SH_mean
dic_SH_anomaly_detrend = signal.detrend(dic_SH_anomaly_detrend)


# Timeseries
fig = plt.figure(figsize = (8, 10))
ax1 = fig.add_axes([0.1,0.7,0.8,0.22])
plt.title('(a) Convection', fontsize = 12)
ax1.set_ylabel('Depth (m)')
ax2.set_ylabel('Degrees C')
plt.xticks([2000,2100, 2200, 2300, 2400, 2500], [' ',' ', ' ', ' ', ' ', ' '])


ax1 = fig.add_axes([0.1,0.42,0.8,0.22])
ax1.plot(heat_pole_total.coord('year').points, (heat_pole_total.data-heat_pole_total_mean)/1e22, color = 'b', lw = 1.5, label = 'Tropics')
ax1.plot(heat_total.coord('year').points, (heat_total.data-heat_mean)/1e22, color = 'k', lw = 1.5, label = 'Global')
plt.title('(b) Heat Content', fontsize = 12)
plt.axhline(0, color = 'k', ls = '--')
plt.ylabel('$10^{22}$ Joules')
plt.ylim([-6, 6])
plt.xticks([2000,2100, 2200, 2300, 2400, 2500], [' ',' ', ' ', ' ', ' ', ' '])


ax1 = fig.add_axes([0.1,0.12,0.8,0.22])
ax1.plot(dic_pole_total.coord('year').points, (dic_pole)/1e15, color = 'b', lw = 1.5, label = 'Southern Ocean')
ax1.plot(dic_pole_total.coord('year').points, (global_carbon)/1e15, color = 'k', lw = 1.5, label = 'Global')
plt.title('(c) Carbon Content', fontsize = 12)
plt.ylabel('Pg Carbon')
plt.axhline(0, color = 'k', ls = '--')
plt.legend(fontsize = 11)
plt.ylim([-5, 5])
ax1.set_xlabel('Time (years)')




plt.show()





