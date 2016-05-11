"""
    Script to plot ohc for different aredi simulations
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


## Calculate OHC ###
names = {'aredi_400', 'aredi_800', 'aredi_2400', 'low_gm'}

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

fig1 = plt.figure(figsize = (6, 8))
ax1 = fig1.add_axes([0.1, 0.79, 0.8, 0.15])
ax1.plot((global_heat_anomaly['aredi_400'])/1e22, color = 'g', label = 'aredi = 400 m$^2$s$^{-1}$')
ax1.plot(heat_SH_anomaly['aredi_400'].data/1e22, color = 'g', ls = '--')

ax2 = fig1.add_axes([0.1, 0.56, 0.8, 0.15])
ax2.plot(global_heat_anomaly['aredi_800']/1e22, color = 'b', label = 'aredi = 800 m$^2$s$^{-1}$')
ax2.plot(heat_SH_anomaly['aredi_800'].data/1e22, color = 'b', ls = '--')

ax3 = fig1.add_axes([0.1, 0.33, 0.8, 0.15])
ax3.plot(global_heat_anomaly['aredi_2400']/1e22, color = 'r', label = 'aredi = 2400 m$^2$s$^{-1}$')
ax3.plot(heat_SH_anomaly['aredi_2400'].data/1e22, color = 'r', ls = '--')
ax3.set_xlim([0,500])

ax4 = fig1.add_axes([0.1, 0.1, 0.8, 0.15])
ax4.plot(global_heat_anomaly['low_gm']/1e22, color = 'k')
ax4.plot(heat_SH_anomaly['low_gm'].data/1e22, color = 'k', ls = '--')
ax4.set_xlim([0,500])

## Detrended
fig1 = plt.figure(figsize = (6, 8))
ax1 = fig1.add_axes([0.15, 0.79, 0.75, 0.15])
ax1.plot((global_heat_detrend['aredi_400'])/1e22, color = 'k', label = 'aredi = 400 m$^2$s$^{-1}$')
ax1.plot(heat_SH_detrend['aredi_400']/1e22, color = 'k', ls = '--')
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11)

ax2 = fig1.add_axes([0.15, 0.56, 0.75, 0.15])
ax2.plot(global_heat_detrend['aredi_800']/1e22, color = 'k', label = 'aredi = 800 m$^2$s$^{-1}$')
ax2.plot(heat_SH_detrend['aredi_800']/1e22, color = 'k', ls = '--')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)
plt.ylabel('10$^{22}$ Joules')

ax3 = fig1.add_axes([0.15, 0.33, 0.75, 0.15])
ax3.plot(global_heat_detrend['aredi_2400']/1e22, color = 'k', label = 'aredi = 2400 m$^2$s$^{-1}$')
ax3.plot(heat_SH_detrend['aredi_2400']/1e22, color = 'k', ls = '--')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])

ax3 = fig1.add_axes([0.15, 0.1, 0.75, 0.15])
ax3.plot(global_heat_detrend['low_gm']/1e22, color = 'k')
ax3.plot(heat_SH_detrend['low_gm']/1e22, color = 'k', ls = '--')
plt.title('(c) GM$_{min}$ = 600 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])
plt.savefig('notes/aredi_ohc_timeseries.png', bbox_inches='tight')


# Double Check the linear detrend method
x = np.arange(0, 500, 1)
x2 = np.arange(0, 507, 1)
p1 = np.poly1d(np.polyfit(x, global_heat_anomaly['aredi_400']/1e22, 1))
p2 = np.poly1d(np.polyfit(x, global_heat_anomaly['aredi_800']/1e22, 1))
p3 = np.poly1d(np.polyfit(x2, global_heat_anomaly['aredi_2400']/1e22, 1))
p7 = np.poly1d(np.polyfit(x2, global_heat_anomaly['low_gm']/1e22, 1))


p4 = np.poly1d(np.polyfit(x, global_heat_anomaly['aredi_400']/1e22, 4))
p5 = np.poly1d(np.polyfit(x, global_heat_anomaly['aredi_800']/1e22, 4))
p6 = np.poly1d(np.polyfit(x2, global_heat_anomaly['aredi_2400']/1e22, 4))
p8 = np.poly1d(np.polyfit(x2, global_heat_anomaly['low_gm']/1e22, 4))


fig1 = plt.figure(figsize = (6, 8))
ax1 = fig1.add_axes([0.1, 0.79, 0.8, 0.15])
ax1.plot((global_heat_anomaly['aredi_400'])/1e22, color = 'k', label = 'aredi = 400 m$^2$s$^{-1}$')
plt.plot(x, p1(x), color = 'k', ls = '--')
plt.plot(x, p4(x), color = 'b', ls = '--')


ax2 = fig1.add_axes([0.1, 0.56, 0.8, 0.15])
ax2.plot(global_heat_anomaly['aredi_800']/1e22, color = 'k', label = 'aredi = 800 m$^2$s$^{-1}$')
plt.plot(x, p2(x), color = 'k', ls = '--')
plt.plot(x, p5(x), color = 'b', ls = '--')


ax3 = fig1.add_axes([0.1, 0.33, 0.8, 0.15])
ax3.plot(global_heat_anomaly['aredi_2400']/1e22, color = 'k', label = 'aredi = 2400 m$^2$s$^{-1}$')
plt.plot(x, p3(x), color = 'k', ls = '--')
plt.plot(x, p6(x), color = 'b', ls = '--')

ax3 = fig1.add_axes([0.1, 0.1, 0.8, 0.15])
ax3.plot(global_heat_anomaly['low_gm']/1e22, color = 'k')
plt.plot(x, p7(x), color = 'k', ls = '--')
plt.plot(x, p8(x), color = 'b', ls = '--')

ax3.set_xlim([0,500])




## Detrend with polynomial
p1 = np.poly1d(np.polyfit(x, global_heat_anomaly['aredi_400']/1e22, 4))
p2 = np.poly1d(np.polyfit(x, global_heat_anomaly['aredi_800']/1e22, 4))
p3 = np.poly1d(np.polyfit(x2, global_heat_anomaly['aredi_2400']/1e22, 4))
p7 = np.poly1d(np.polyfit(x2, global_heat_anomaly['low_gm']/1e22, 4))


p4 = np.poly1d(np.polyfit(x, heat_SH_anomaly['aredi_400'].data/1e22, 4))
p5 = np.poly1d(np.polyfit(x, heat_SH_anomaly['aredi_800'].data/1e22, 4))
p6 = np.poly1d(np.polyfit(x2, heat_SH_anomaly['aredi_2400'].data/1e22, 4))
p8 = np.poly1d(np.polyfit(x2, heat_SH_anomaly['low_gm'].data/1e22, 4))


y1 = global_heat_anomaly['aredi_400']/1e22 - p1(x)
y2 = global_heat_anomaly['aredi_800']/1e22 - p2(x)
y3 = global_heat_anomaly['aredi_2400']/1e22 - p3(x2)
y7 = global_heat_anomaly['low_gm']/1e22 - p7(x2)


y4 = heat_SH_anomaly['aredi_400'].data/1e22 - p4(x)
y5 = heat_SH_anomaly['aredi_800'].data/1e22 - p5(x)
y6 = heat_SH_anomaly['aredi_2400'].data/1e22 - p6(x2)
y8 = heat_SH_anomaly['low_gm'].data/1e22 - p8(x2)


fig1 = plt.figure(figsize = (6, 8))
ax1 = fig1.add_axes([0.15, 0.79, 0.75, 0.15])
ax1.plot(y1, color = 'k', label = 'aredi = 400 m$^2$s$^{-1}$')
plt.plot(y4, color = 'k', ls = '--')
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11)
plt.ylim([-6, 6])

ax2 = fig1.add_axes([0.15, 0.56, 0.75, 0.15])
ax2.plot(y2, color = 'k', label = 'aredi = 800 m$^2$s$^{-1}$')
plt.plot(y5, color = 'k', ls = '--')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)
plt.ylim([-6, 6])

ax3 = fig1.add_axes([0.15, 0.33, 0.75, 0.15])
ax3.plot(y3, color = 'k', label = 'aredi = 2400 m$^2$s$^{-1}$')
plt.plot(y6, color = 'k', ls = '--')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])
plt.ylabel('                                 10$^{22}$ Joules')
plt.ylim([-6, 6])


ax3 = fig1.add_axes([0.15, 0.1, 0.75, 0.15])
ax3.plot(y7, color = 'k')
plt.plot(y8, color = 'k', ls = '--')
plt.title('(c) GM$_{min}$ = 600 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])
plt.xlabel('Time [years]')


plt.savefig('notes/aredi_ohc_timeseries_quadratic_detrend.png', bbox_inches='tight')
plt.ylim([-6, 6])

plt.show()







