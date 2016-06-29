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
from scipy.optimize import curve_fit
from scipy import log as log
from scipy import exp as exp

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
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11)

ax2 = fig1.add_axes([0.1, 0.56, 0.8, 0.15])
ax2.plot(global_heat_anomaly['aredi_800']/1e22, color = 'b', label = 'aredi = 800 m$^2$s$^{-1}$')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)

ax3 = fig1.add_axes([0.1, 0.33, 0.8, 0.15])
ax3.plot(global_heat_anomaly['aredi_2400']/1e22, color = 'r', label = 'aredi = 2400 m$^2$s$^{-1}$')
ax3.set_xlim([0,500])
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
plt.ylabel('                                 10$^{22}$ Joules')


ax4 = fig1.add_axes([0.1, 0.1, 0.8, 0.15])
ax4.plot(global_heat_anomaly['low_gm']/1e22, color = 'k')
ax4.set_xlim([0,500])
plt.title('(c) GM$_{min}$ = 600 m$^2$s$^{-1}$', fontsize = 11)
plt.xlabel('Time [years]')



# Fit logarithmic function

yarr =  global_heat_anomaly['low_gm'][-400:].data/1e22
xarr = np.arange(0, 400, 1)
xarr2 = np.arange(0, 1000, 1)


def func(x, a, b, c, d):
    
    return a*x**3 + b*x**2 + c*x +d

popt, pcov = curve_fit(func, xarr, yarr)



plt.figure()
plt.plot(xarr, yarr, '*')
plt.plot(xarr2, func(xarr2, *popt), 'g')
plt.ylim([-60,100])


plt.show()







