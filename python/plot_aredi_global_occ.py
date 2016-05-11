"""
    Script to plot occ for different aredi simulations
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

for name in names:
    carbon[name], carbon_sumz[name], carbon_global[name] = calc_occ(dic[name], rhodzt[name], area[name])
    global_carbon_mean[name] = carbon_global[name].collapsed('time', iris.analysis.MEAN)
    global_carbon_anomaly[name] = carbon_global[name].data-global_carbon_mean[name].data
    global_carbon_detrend[name] = signal.detrend(global_carbon_anomaly[name])

    # Southern Hemisphere:
    constraint = iris.Constraint(latitude=lambda y: -90 < y < -50)
    carbon_SH_tmp = carbon_sumz[name].extract(constraint)
    area_SH = area[name].extract(constraint)
    b = carbon_SH_tmp*area_SH
    carbon_SH[name] = b.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    mean_SH = carbon_SH[name].collapsed('time', iris.analysis.MEAN)
    carbon_SH_anomaly[name] = carbon_SH[name] - mean_SH
    carbon_SH_detrend[name] = signal.detrend(carbon_SH_anomaly[name].data)


fig1 = plt.figure(figsize = (6, 8))
ax1 = fig1.add_axes([0.1, 0.79, 0.8, 0.15])
ax1.plot((global_carbon_anomaly['aredi_400'])/1e15, color = 'g', label = 'aredi = 400 m$^2$s$^{-1}$')
ax1.plot(carbon_SH_anomaly['aredi_400'].data/1e15, color = 'g', ls = '--')

ax2 = fig1.add_axes([0.1, 0.56, 0.8, 0.15])
ax2.plot(global_carbon_anomaly['aredi_800']/1e15, color = 'b', label = 'aredi = 800 m$^2$s$^{-1}$')
ax2.plot(carbon_SH_anomaly['aredi_800'].data/1e15, color = 'b', ls = '--')

ax3 = fig1.add_axes([0.1, 0.33, 0.8, 0.15])
ax3.plot(global_carbon_anomaly['aredi_2400']/1e15, color = 'r', label = 'aredi = 2400 m$^2$s$^{-1}$')
ax3.plot(carbon_SH_anomaly['aredi_2400'].data/1e15, color = 'r', ls = '--')
ax3.set_xlim([0,500])

ax3 = fig1.add_axes([0.1, 0.1, 0.8, 0.15])
ax3.plot(global_carbon_anomaly['low_gm']/1e15, color = 'k')
ax3.plot(carbon_SH_anomaly['low_gm'].data/1e15, color = 'k', ls = '--')
ax3.set_xlim([0,500])

## Detrended
fig1 = plt.figure(figsize = (6, 8))
ax1 = fig1.add_axes([0.15, 0.79, 0.75, 0.15])
ax1.plot((global_carbon_detrend['aredi_400'])/1e15, color = 'k', label = 'aredi = 400 m$^2$s$^{-1}$')
ax1.plot(carbon_SH_detrend['aredi_400']/1e15, color = 'k', ls = '--')
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11)

ax2 = fig1.add_axes([0.15, 0.56, 0.75, 0.15])
ax2.plot(global_carbon_detrend['aredi_800']/1e15, color = 'k', label = 'aredi = 800 m$^2$s$^{-1}$')
ax2.plot(carbon_SH_detrend['aredi_800']/1e15, color = 'k', ls = '--')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)
plt.ylabel('PgC')

ax3 = fig1.add_axes([0.15, 0.33, 0.75, 0.15])
ax3.plot(global_carbon_detrend['aredi_2400']/1e15, color = 'k', label = 'aredi = 2400 m$^2$s$^{-1}$')
ax3.plot(carbon_SH_detrend['aredi_2400']/1e15, color = 'k', ls = '--')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])

ax3 = fig1.add_axes([0.15, 0.1, 0.75, 0.15])
ax3.plot(global_carbon_detrend['low_gm']/1e15, color = 'k')
ax3.plot(carbon_SH_detrend['low_gm']/1e15, color = 'k', ls = '--')
plt.title('(c) GM$_{min}$ = 600 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])
plt.savefig('notes/aredi_occ_timeseries.png', bbox_inches='tight')


# Double Check the linear detrend method
x = np.arange(0, 500, 1)
x2 = np.arange(0, 507, 1)
p1 = np.polyfit(x, global_carbon_anomaly['aredi_400']/1e15, 1)
p2 = np.polyfit(x, global_carbon_anomaly['aredi_800']/1e15, 1)
p3 = np.polyfit(x2, global_carbon_anomaly['aredi_2400']/1e15, 1)
p7 = np.polyfit(x2, global_carbon_anomaly['low_gm']/1e15, 1)
p4 = np.poly1d(np.polyfit(x2, global_carbon_anomaly['aredi_2400']/1e15, 4))
p5 = np.poly1d(np.polyfit(x, global_carbon_anomaly['aredi_400']/1e15, 4))
p6 = np.poly1d(np.polyfit(x, global_carbon_anomaly['aredi_800']/1e15, 4))
p8 = np.poly1d(np.polyfit(x2, global_carbon_anomaly['low_gm']/1e15, 4))



fig1 = plt.figure(figsize = (6, 8))
ax1 = fig1.add_axes([0.15, 0.79, 0.75, 0.15])
ax1.plot((global_carbon_anomaly['aredi_400'])/1e15, color = 'k', label = 'aredi = 400 m$^2$s$^{-1}$')
plt.plot(x, p1[0]*x + p1[1], color = 'k', ls = '--')
plt.plot(x, p5(x), color = 'b', ls = '--')
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11)

ax2 = fig1.add_axes([0.15, 0.56, 0.75, 0.15])
ax2.plot(global_carbon_anomaly['aredi_800']/1e15, color = 'k', label = 'aredi = 800 m$^2$s$^{-1}$')
plt.plot(x, p2[0]*x + p2[1], color = 'k', ls = '--')
plt.plot(x, p6(x), color = 'b', ls = '--')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)
plt.ylabel('PgC')

ax3 = fig1.add_axes([0.15, 0.33, 0.75, 0.15])
ax3.plot(global_carbon_anomaly['aredi_2400']/1e15, color = 'k', label = 'aredi = 2400 m$^2$s$^{-1}$')
plt.plot(x, p3[0]*x + p3[1], color = 'k', ls = '--')
plt.plot(x2, p4(x2), color = 'b', ls = '--')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])

ax3 = fig1.add_axes([0.15, 0.1, 0.75, 0.15])
ax3.plot(global_carbon_anomaly['low_gm']/1e15, color = 'k')
plt.plot(x, p7[0]*x + p7[1], color = 'k', ls = '--')
plt.plot(x2, p8(x2), color = 'b', ls = '--')
plt.title('(c) GM$_{min}$ = 600 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])


## Detrend with polynomial
p1 = np.poly1d(np.polyfit(x, global_carbon_anomaly['aredi_400']/1e15, 4))
p2 = np.poly1d(np.polyfit(x, global_carbon_anomaly['aredi_800']/1e15, 4))
p3 = np.poly1d(np.polyfit(x2, global_carbon_anomaly['aredi_2400']/1e15, 4))
p7 = np.poly1d(np.polyfit(x2, global_carbon_anomaly['low_gm']/1e15, 4))


p4 = np.poly1d(np.polyfit(x, carbon_SH_anomaly['aredi_400'].data/1e15, 4))
p5 = np.poly1d(np.polyfit(x, carbon_SH_anomaly['aredi_800'].data/1e15, 4))
p6 = np.poly1d(np.polyfit(x2, carbon_SH_anomaly['aredi_2400'].data/1e15, 4))
p8 = np.poly1d(np.polyfit(x2, carbon_SH_anomaly['low_gm'].data/1e15, 4))


y1 = global_carbon_anomaly['aredi_400']/1e15 - p1(x)
y2 = global_carbon_anomaly['aredi_800']/1e15 - p2(x)
y3 = global_carbon_anomaly['aredi_2400']/1e15 - p3(x2)
y7 = global_carbon_anomaly['low_gm']/1e15 - p7(x2)


y4 = carbon_SH_anomaly['aredi_400'].data/1e15 - p4(x)
y5 = carbon_SH_anomaly['aredi_800'].data/1e15 - p5(x)
y6 = carbon_SH_anomaly['aredi_2400'].data/1e15 - p6(x2)
y8 = carbon_SH_anomaly['low_gm'].data/1e15 - p8(x2)


fig1 = plt.figure(figsize = (6, 8))
ax1 = fig1.add_axes([0.15, 0.79, 0.75, 0.15])
ax1.plot(y1, color = 'k', label = 'aredi = 400 m$^2$s$^{-1}$')
plt.plot(y4, color = 'k', ls = '--')
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11)
plt.ylim([-5, 5])

ax2 = fig1.add_axes([0.15, 0.56, 0.75, 0.15])
ax2.plot(y2, color = 'k', label = 'aredi = 800 m$^2$s$^{-1}$')
plt.plot(y5, color = 'k', ls = '--')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)
plt.ylim([-5, 5])

ax3 = fig1.add_axes([0.15, 0.33, 0.75, 0.15])
ax3.plot(y3, color = 'k', label = 'aredi = 2400 m$^2$s$^{-1}$')
plt.plot(y6, color = 'k', ls = '--')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])
plt.ylim([-5, 5])
plt.ylabel('                              PgC')


ax3 = fig1.add_axes([0.15, 0.1, 0.75, 0.15])
ax3.plot(y7, color = 'k')
plt.plot(y8, color = 'k', ls = '--')
plt.title('(c) GM$_{min}$ = 600 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])
plt.ylim([-5, 5])
plt.xlabel('Time [years]')
plt.savefig('notes/aredi_occ_timeseries_quadratic_detrend.png', bbox_inches='tight')


plt.show()







