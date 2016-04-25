###
# Script to calculate and plot ocean carbon content for six different regions in the ocean
###


import iris.analysis.cartography
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from scipy import signal

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_occ
from calc_region import calc_six_regions
from smooth import smooth

# Load Temp Data
PATHS={'aredi_800':'~/control_aredi800/data/newCO2_control_800/'}

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

names = {'aredi_800'}

carbon = {}
carbon_sumz = {}
carbon_global = {}
global_carbon_anomaly = {}
global_carbon_mean = {}

for name in names:
    carbon[name], carbon_sumz[name], carbon_global[name] = calc_occ(dic[name], rhodzt[name], area[name])
    global_carbon_mean[name] = carbon_global[name].collapsed('time', iris.analysis.MEAN)
    global_carbon_anomaly[name] = carbon_global[name].data-global_carbon_mean[name].data

SO, SPac, SAtl, NPac, NAtl, Arctic = calc_six_regions(carbon_sumz['aredi_800'])
a1, a2, a3, a4, a5, a6 = calc_six_regions(area['aredi_800'])

# Sum in each region:
SO = SO*a1
SPac = SPac*a2
SAtl = SAtl*a3
NPac = NPac*a4
NAtl = NAtl*a5
Arctic = Arctic*a6

SO = SO.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
SPac = SPac.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
SAtl = SAtl.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
NPac = NPac.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
NAtl = NAtl.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
Arctic = Arctic.collapsed(['longitude', 'latitude'], iris.analysis.SUM)

SO_mean = SO.collapsed('time', iris.analysis.MEAN)
SPac_mean = SPac.collapsed('time', iris.analysis.MEAN)
SAtl_mean = SAtl.collapsed('time', iris.analysis.MEAN)
NPac_mean = NPac.collapsed('time', iris.analysis.MEAN)
NAtl_mean = NAtl.collapsed('time', iris.analysis.MEAN)
Arctic_mean = Arctic.collapsed('time', iris.analysis.MEAN)

SO = SO-SO_mean
SPac = SPac - SPac_mean
SAtl = SAtl - SAtl_mean
NPac = NPac - NPac_mean
NAtl = NAtl - NAtl_mean
Arctic = Arctic - Arctic_mean

fig = plt.figure()
plt.plot(SO.data/1e15, label = 'Southern Ocean')
plt.plot(SPac.data/1e15, label = 'South Pacific')
plt.plot(SAtl.data/1e15, label = 'South Atlantic')
plt.plot(NPac.data/1e15, label = 'North Pacific')
plt.plot(NAtl.data/1e15, label = 'North Atlantic')
plt.plot(Arctic.data/1e15, label = 'Arctic')
plt.plot(global_carbon_anomaly['aredi_800'].data/1e15, color = 'k', lw = 2)
plt.legend()


# Detrend
SO_detrend = signal.detrend(SO.data)
SPac_detrend = signal.detrend(SPac.data)
SAtl_detrend = signal.detrend(SAtl.data)
NPac_detrend = signal.detrend(NPac.data)
NAtl_detrend = signal.detrend(NAtl.data)
Arctic_detrend = signal.detrend(Arctic.data) 
global_carbon_detrend = signal.detrend(global_carbon_anomaly['aredi_800'].data)

fig = plt.figure()
plt.plot(SO_detrend/1e15, label = 'Southern Ocean')
plt.plot(SPac_detrend/1e15, label = 'South Pacific')
plt.plot(SAtl_detrend/1e15, label = 'South Atlantic')
plt.plot(NPac_detrend/1e15, label = 'North Pacific')
plt.plot(NAtl_detrend/1e15, label = 'North Atlantic')
plt.plot(Arctic_detrend/1e15, label = 'Arctic')
plt.plot(global_carbon_detrend/1e15, color = 'k', lw = 2)
plt.legend()
plt.xlabel('Time [years]')
plt.ylabel('Carbon content anomaly [Pg C]')

x = 21
SO_smooth = smooth(SO_detrend, window_len=x, window='hanning')
SPac_smooth = smooth(SPac_detrend, window_len=x, window='hanning')
SAtl_smooth = smooth(SAtl_detrend, window_len=x, window='hanning')
NPac_smooth = smooth(NPac_detrend, window_len=x, window='hanning')
NAtl_smooth = smooth(NAtl_detrend, window_len=x, window='hanning')
Arctic_smooth = smooth(Arctic_detrend, window_len=x, window='hanning')
global_carbon_smooth = smooth(global_carbon_detrend, window_len = x, window='hanning')

fig = plt.figure()
plt.plot(SO_smooth/1e15, label = 'Southern Ocean')
plt.plot(SPac_smooth/1e15, label = 'South Pacific')
plt.plot(SAtl_smooth/1e15, label = 'South Atlantic')
plt.plot(NPac_smooth/1e15, label = 'North Pacific')
plt.plot(NAtl_smooth/1e15, label = 'North Atlantic')
plt.plot(Arctic_smooth/1e15, label = 'Arctic')
plt.plot(global_carbon_smooth/1e15, color = 'k', lw = 2)

plt.text(390,3.6, 'Southern Ocean', color = 'b')
plt.text(390,3.3, 'South Pacific', color = 'g')
plt.text(390,3.0, 'South Atlantic', color = 'r')
plt.text(390,2.7, 'North Pacific', color = 'c')
plt.text(390,2.4, 'North Atlantic', color = 'purple')
plt.text(390,2.1, 'Arctic', color = 'gold')
plt.text(390,1.8, 'Global', color = 'k')
plt.xlim([0, 500])
plt.xlabel('Time [years]')
plt.ylabel('Carbon content anomaly [Pg C]')

fig = plt.figure()
sum = (SO_smooth + SPac_smooth + SAtl_smooth + NPac_smooth + NAtl_smooth + Arctic_smooth)/1e15
plt.plot(SO_smooth/1e15, label = 'Southern Ocean')
plt.plot(SPac_smooth/1e15, label = 'South Pacific')
plt.plot(SAtl_smooth/1e15, label = 'South Atlantic')
plt.plot(NPac_smooth/1e15, label = 'North Pacific')
plt.plot(NAtl_smooth/1e15, label = 'North Atlantic')
plt.plot(Arctic_smooth/1e15, label = 'Arctic')
plt.plot(global_carbon_smooth/1e15, color = 'k', lw = 2)
plt.plot(sum, color = 'k', ls = '--', lw = 2)
plt.xlim([100, 200])
plt.xlabel('Time [years]')
plt.ylabel('Carbon content anomaly [Pg C]')
plt.show()
