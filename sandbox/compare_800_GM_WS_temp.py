
"""
    Script to compare the 1500-2000m Weddell Sea Temperature for Aredi=800 and GM_min=600 simulations
"""

import iris
import iris.analysis.cartography
import matplotlib.pyplot as plt
import numpy as np
import glob
import time
import numpy as np
import matplotlib.mlab as mlab
from scipy.stats import normaltest
from scipy.stats import genextreme

import sys # access system routines
sys.path.append('/home/jthom143/python_functions')
import colormaps as cmaps


import sys # access system routines
sys.path.append('~/control_aredi800/python')
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

temp = {}
mld = {}
mld_ws = {}
temp_ws = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    temp[name] = iris.load_cube(PATH+'temp.nc')
    mld[name] = iris.load_cube(PATH+'mld.nc')
    
    temp_ws[name]= calc_weddell_sea(temp[name])
    mld_ws[name] = calc_weddell_sea(mld[name])
    
    if name == 'aredi_800':
        temp[name] = temp[name][-500:]
    if name == 'aredi_2400':
        temp[name] = temp[name][:500]
        mld[name] = mld[name][:500]
    if name == 'low_gm':
        temp[name] = temp[name][:500]
        mld[name] = mld[name][:500]

    mld_ws[name].coord('latitude').guess_bounds()
    mld_ws[name].coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mld_ws[name])
    mld_ws[name] = mld_ws[name].collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    
    temp_ws[name].coord('latitude').guess_bounds()
    temp_ws[name].coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(temp_ws[name])
    temp_ws[name] = temp_ws[name].collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)

# Detrend the low_gm case
x2 = np.arange(0, 507, 1)
p = np.poly1d(np.polyfit(x2[50:], temp_ws['low_gm'][50:,18].data, 4))
y = temp_ws['low_gm'][50:,18].data - p(x2[50:]) + np.mean(temp_ws['low_gm'][50:,18].data)



# Create Figure
plt.figure()
time = np.arange(0, 500, 1)
time2 = np.arange(0, 507, 1)
plt.plot(time, temp_ws['aredi_400'][:,18].data, color = 'b', lw = 1.5, label = 'Aredi = 400')
plt.plot(time, temp_ws['aredi_800'][:,18].data, color = 'g', lw = 1.5, label = 'Aredi = 800')
plt.plot(time2, temp_ws['aredi_2400'][:,18].data, color = 'r', lw = 1.5, label = 'Aredi = 2400')
plt.plot(time2[50:], y, color = 'g', lw = 1.5, ls = '--', label = 'GM$_{min}$ = 600')
plt.xlim([0,500])
plt.xlabel('Time [years]')
plt.ylabel('Temperature [$^o$C]')
plt.title('1400m Weddell Sea Temperature')
plt.legend(fontsize = 10)
plt.axhline(0.76327787689, color = 'k', ls = '--', lw = 1.5)


# PDFS


# Normal test
test = normaltest(temp_ws['aredi_2400'][:,18].data)

plt.figure()
plt.subplot(2,2,1)
plt.hist(temp_ws['aredi_400'][:,18].data, 50, color = 'b', normed=True)
plt.ylabel('Probability Density')
plt.title('Aredi = 400 m$^2$s$^{-1}$', fontsize = 11)
x = np.linspace(0,2.5,100)

param = genextreme.fit(temp_ws['aredi_400'][:,18].data) # distribution fitting
pdf_fitted = genextreme.pdf(x,c = param[0], loc=param[1],scale=param[2])
plt.plot(x,pdf_fitted,'k--', lw = 1.5)

plt.subplot(2,2,2)
plt.hist(temp_ws['aredi_800'][:,18].data, 50, color = 'g', normed=True)
plt.title('Aredi = 800 m$^2$s$^{-1}$', fontsize = 11)

param = genextreme.fit(temp_ws['aredi_800'][:,18].data) # distribution fitting
pdf_fitted_800 = genextreme.pdf(x,c = param[0], loc=param[1],scale=param[2])
plt.plot(x,pdf_fitted_800,'k--', lw = 1.5)

plt.subplot(2,2,3)
n, bins_2400, patches = plt.hist(temp_ws['aredi_2400'][:,18].data, 25, normed=True, color = 'r', label = 'Aredi = 2400 m$^2$s$^{-1}$')
mu = np.mean(temp_ws['aredi_2400'][:,18].data)
sigma = np.std(temp_ws['aredi_2400'][:,18].data)
y_hist_2400 = mlab.normpdf( bins_2400, mu, sigma)
l = plt.plot(bins_2400, y_hist_2400, 'k--', linewidth=1.5)
plt.xlabel('Temperature [$^o$C]')
plt.ylabel('Probability Density')
plt.title('Aredi = 2400 m$^2$s$^{-1}$', fontsize = 11)


plt.subplot(2,2,4)
n, bins, patches = plt.hist(y, 15, color = 'c', normed=True)
plt.xlabel('Temperature [$^o$C]')
plt.title('GM$_{min}$ = 600 m$^2$s$^{-1}$', fontsize = 11)
plt.xlim([2.0,2.5])

mu = np.mean(y)
sigma = np.std(y)
y_hist_gm = mlab.normpdf( x, mu, sigma)
l = plt.plot(x, y_hist_gm, 'k--', linewidth=1.5)

plt.figure()
plt.plot(x,pdf_fitted,'b-', lw = 1.5)
plt.plot(x,pdf_fitted_800,'g-', lw = 1.5)
plt.ylabel('Probability Density')


plt.plot(bins_2400, y_hist_2400, 'r', linewidth=1.5)
l = plt.plot(x, y_hist_gm, 'c-', linewidth=1.5)
plt.axvline(0.76327787689, color = 'k', ls = '--', lw = 1.5)


plt.xlim([0,2.5])
plt.xlabel('Temperature [$^o$C]')


plt.show()

