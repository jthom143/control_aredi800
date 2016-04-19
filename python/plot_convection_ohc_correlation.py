"""
    Script to plot the correlation between convection and heat content
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
#from scipy.stats import genextreme as dist
from scipy.stats import gamma as dist

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_ohc
from calc_region import calc_weddell_sea

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


## Calculate mld in Weddell Sea ##
mld_ws = {}

for name in names:
    mld_ws[name] = calc_weddell_sea(mld[name])
    mld_ws[name].coord('latitude').guess_bounds()
    mld_ws[name].coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mld_ws[name])
    mld_ws[name] = mld_ws[name].collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=grid_areas)

# Plot Weddell Sea MLD
## Detrended
fig1 = plt.figure(figsize = (6, 8))
ax1 = fig1.add_axes([0.15, 0.75, 0.75, 0.15])
ax1.plot(mld_ws['aredi_400'].data, color = 'k', label = 'aredi = 400 m$^2$s$^{-1}$')
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11)

ax2 = fig1.add_axes([0.15, 0.5, 0.75, 0.15])
ax2.plot(mld_ws['aredi_800'].data, color = 'k', label = 'aredi = 800 m$^2$s$^{-1}$')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)
plt.ylabel('Weddell Sea MLD [m]')

ax3 = fig1.add_axes([0.15, 0.25, 0.75, 0.15])
ax3.plot(mld_ws['aredi_2400'].data, color = 'k', label = 'aredi = 2400 m$^2$s$^{-1}$')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
ax3.set_xlim([0,500])
plt.savefig('notes/figures/aredi_mld_timeseries.png', bbox_inches='tight')



'''
fig2 = plt.figure()
data1 = mld_ws['aredi_400'].data
x = np.arange(0, 810, 10)
p1 = dist.fit(data1)
plt.hist(data1, x, normed=True)
plt.plot(x, dist.pdf(x, *p1), lw = 2, color = 'k')
plt.title('Aredi = 400')
plt.xlabel('Mixed Layer Depth [m]')
plt.ylabel('Density')

fig3 = plt.figure()
data2 = mld_ws['aredi_800'].data
p2 = dist.fit(data2)
plt.hist(data2, x, normed=True)
plt.plot(x, dist.pdf(x, *p2), lw = 2, color = 'k')
plt.title('Aredi = 800')
plt.xlabel('Mixed Layer Depth [m]')
plt.ylabel('Density')

fig4 = plt.figure()
data3 = mld_ws['aredi_2400'].data
p3 = dist.fit(data3)
x2 = np.arange(0, 510, 10)
plt.hist(data3, x2, normed=True)
plt.plot(x2, dist.pdf(x2, *p3), lw = 2, color = 'k')
plt.title('Aredi = 2400')
plt.xlabel('Mixed Layer Depth [m]')
plt.ylabel('Density')



fig5 = plt.figure()
plt.plot(x, dist.pdf(x, *p1), lw = 2, color = 'g', label = 'Aredi = 400')
plt.plot(x, dist.pdf(x, *p2), lw = 2, color = 'b', label = 'Aredi = 800')
plt.plot(x2, dist.pdf(x2, *p3), lw = 2, color = 'r', label = 'Aredi = 2400')
plt.xlabel('Mixed Layer Depth [m]')
plt.ylabel('Density')
plt.legend()
plt.title('Gamma Distribution')
plt.savefig('notes/figures/aredi_mld_GammaDist.png', bbox_inches='tight')

fig6 = plt.figure()
plt.plot(x, dist.cdf(x, *p1), lw = 2, color = 'g', label = 'Aredi = 400')
plt.plot(x, dist.cdf(x, *p2), lw = 2, color = 'b', label = 'Aredi = 800')
plt.plot(x2, dist.cdf(x2, *p3), lw = 2, color = 'r', label = 'Aredi = 2400')
plt.xlabel('Mixed Layer Depth [m]')
plt.ylabel('Density')
plt.legend()
plt.title('Gamma Cumulative Distribution')
#plt.savefig('notes/figures/aredi_mld_GammaDist.png', bbox_inches='tight')
'''

### Correlate Heat Content with MLD
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

# Save variables for later use
for name in names:
    np.save(name+'mld', mld_ws[name].data)
    np.save(name+'global_heat', global_heat_detrend[name])
    np.save(name+'sh_heat', heat_SH_detrend[name])

'''
# Find the maximum correlation
x = {}
y = {}
y2 = {}
xcorr = {}
min_corr = {}
min_lag = {}

for name in names:
    x[name] = mld_ws[name].data
    y[name] = global_heat_detrend[name]/1e22
    y2[name] = heat_SH_detrend[name]/1e22

    xcorr[name] = np.correlate(x[name], y[name], mode="full")
    a = np.arange(xcorr[name].size)
    lags = a - (x[name].size-1)

    min_corr[name] = np.argmin(xcorr[name])
    min_lag[name] = lags[min_corr[name]]

fig9 = plt.figure()
plt.xcorr(x['aredi_400'], y['aredi_400'], maxlags = 100, usevlines = False)
plt.xcorr(x['aredi_400'], y2['aredi_400'], maxlags = 100, usevlines = False)

fig9 = plt.figure()
plt.xcorr(x['aredi_800'], y['aredi_800'], maxlags = 100, usevlines = False)
plt.xcorr(x['aredi_800'], y2['aredi_800'], maxlags = 100, usevlines = False)

fig9 = plt.figure()
plt.xcorr(x['aredi_2400'], y['aredi_2400'], maxlags = 100, usevlines = False)
plt.xcorr(x['aredi_2400'], y2['aredi_2400'], maxlags = 100, usevlines = False)
'''
plt.show()
