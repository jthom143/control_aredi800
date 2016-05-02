###
# Script to calculate and plot ocean heat content for six different regions in the ocean
###


import iris.analysis.cartography
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from scipy import signal

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_ohc
from calc_region import calc_six_regions, calc_zonal_regions
from smooth import smooth

# Load Temp Data
PATHS={'aredi_800':'~/control_aredi800/data/newCO2_control_800/'}

temp = {}
rhodzt = {}
area = {}


print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    temp[name] = iris.load_cube(PATH+'temp.nc')
    rhodzt[name] = iris.load_cube(PATH+'rho_dzt.nc')
    area[name] = iris.load_cube(PATH+'area_t.nc')


names = ['aredi_800']

heat = {}
heat_sumz = {}
heat_global = {}
global_heat_anomaly = {}
global_heat_mean = {}

for name in names:
    heat[name], heat_sumz[name], heat_global[name] = calc_ohc(temp[name], rhodzt[name], area[name])
    global_heat_mean[name] = heat_global[name].collapsed('time', iris.analysis.MEAN)
    global_heat_anomaly[name] = heat_global[name].data-global_heat_mean[name].data

#SO, SPac, SAtl, NPac, NAtl, Arctic = calc_six_regions(heat_sumz['aredi_800'])
#a1, a2, a3, a4, a5, a6 = calc_six_regions(area['aredi_800'])

SO, SMid, NMid, arctic, tropics =  calc_zonal_regions(heat_sumz['aredi_800'])
a1, a2, a3, a4, a5 = calc_zonal_regions(area['aredi_800'])

# Sum in each region:
SO = SO*a1
SMid = SMid*a2
NMid = NMid*a3
arctic = arctic*a4
tropics = tropics*a5

SO = SO.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
SMid = SMid.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
NMid = NMid.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
tropics = tropics.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
arctic = arctic.collapsed(['longitude', 'latitude'], iris.analysis.SUM)

SO_mean = SO.collapsed('time', iris.analysis.MEAN)
SMid_mean = SMid.collapsed('time', iris.analysis.MEAN)
NMid_mean = NMid.collapsed('time', iris.analysis.MEAN)
tropics_mean = tropics.collapsed('time', iris.analysis.MEAN)
arctic_mean = arctic.collapsed('time', iris.analysis.MEAN)

SO = SO-SO_mean
SMid = SMid - SMid_mean
NMid = NMid - NMid_mean
tropics = tropics - tropics_mean
arctic = arctic - arctic_mean

fig = plt.figure()
plt.plot(SO.data/1e22, label = 'Southern Ocean')
plt.plot(SMid.data/1e22, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid.data/1e22, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics.data/1e22, label = 'Tropics')
plt.plot(arctic.data/1e22, label = 'Arctic')
plt.plot(global_heat_anomaly['aredi_800'].data/1e22, color = 'k', lw = 2)
plt.legend()


# Detrend
SO_detrend = signal.detrend(SO.data)
SMid_detrend = signal.detrend(SMid.data)
NMid_detrend = signal.detrend(NMid.data)
tropics_detrend = signal.detrend(tropics.data)
arctic_detrend = signal.detrend(arctic.data)
global_heat_detrend = signal.detrend(global_heat_anomaly['aredi_800'].data)

fig = plt.figure()
plt.plot(SO_detrend/1e22, label = 'Southern Ocean')
plt.plot(SMid_detrend/1e22, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid_detrend/1e22, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics_detrend/1e22, label = 'Tropics')
plt.plot(arctic_detrend/1e22, label = 'Arctic')
plt.plot(global_heat_detrend/1e22, color = 'k', lw = 2)
plt.legend()
plt.xlabel('Time [years]')
plt.ylabel('Heat content anomaly [10$^{22}$ J]')

x = 21
SO_smooth = smooth(SO_detrend, window_len=x, window='hanning')
SMid_smooth = smooth(SMid_detrend, window_len=x, window='hanning')
NMid_smooth = smooth(NMid_detrend, window_len=x, window='hanning')
tropics_smooth = smooth(tropics_detrend, window_len=x, window='hanning')
arctic_smooth = smooth(arctic_detrend, window_len=x, window='hanning')
global_heat_smooth = smooth(global_heat_detrend, window_len = x, window='hanning')

fig = plt.figure(figsize=(12,8))
plt.plot(SO_smooth/1e22, label = 'Southern Ocean')
plt.plot(SMid_smooth/1e22, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid_smooth/1e22, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics_smooth/1e22, label = 'Tropics')
plt.plot(arctic_smooth/1e22, label = 'Arctic')
plt.plot(global_heat_smooth/1e22, color = 'k', lw = 2)

plt.text(390,5.6, 'Southern Ocean', color = 'b')
plt.text(390,5.3, 'SH Mid-Lats', color = 'g')
plt.text(390,5.0, 'NH Mid-Lats', color = 'r')
plt.text(390,5.7, 'Tropics', color = 'c')
plt.text(390,5.4, 'Arctic', color = 'purple')
plt.text(390,5.1, 'Global', color = 'k')
plt.xlim([0, 500])
plt.xlabel('Time [years]')
plt.ylabel('Heat content anomaly [10$^{22}$ J]')
plt.savefig('/home/jthom143/control_aredi800/notes/figures/ohc_regions.png')

fig = plt.figure()
sum = (SO_smooth + SMid_smooth +  NMid_smooth + tropics_smooth + arctic_smooth)/1e22
plt.plot(SO_smooth/1e22, label = 'Southern Ocean')
plt.plot(SMid_smooth/1e22, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid_smooth/1e22, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics_smooth/1e22, label = 'Tropics')
plt.plot(arctic_smooth/1e22, label = 'Arctic')

plt.plot(global_heat_smooth/1e22, color = 'k', lw = 2)
plt.plot(sum, color = 'k', ls = '--', lw = 2)
plt.xlim([100, 200])
plt.xlabel('Time [years]')
plt.ylabel('Heat content anomaly [10$^{22}$ J]')
plt.show()

