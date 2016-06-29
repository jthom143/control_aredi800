###
# Script to calculate and plot ocean carbon and heat content for six different regions in the ocean
###


import iris.analysis.cartography
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from scipy import signal

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_occ, calc_ohc
from calc_region import calc_six_regions, calc_zonal_regions
from smooth import smooth

# Load Temp Data
PATHS={'aredi_800':'~/control_aredi800/data/newCO2_control_800/'}

dic = {}
temp = {}
rhodzt = {}
area = {}


print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    dic[name] = iris.load_cube(PATH+'dic.nc')
    temp[name] = iris.load_cube(PATH+'temp.nc')
    rhodzt[name] = iris.load_cube(PATH+'rho_dzt.nc')
    area[name] = iris.load_cube(PATH+'area_t.nc')
    
    if name == 'aredi_800':
        dic[name] = dic[name][-500:]

names = ['aredi_800']

carbon = {}
carbon_sumz = {}
carbon_global = {}
global_carbon_anomaly = {}
global_carbon_mean = {}

heat = {}
heat_sumz = {}
heat_global = {}
global_heat_anomaly = {}
global_heat_mean = {}

for name in names:
    carbon[name], carbon_sumz[name], carbon_global[name] = calc_occ(dic[name], rhodzt[name], area[name])
    global_carbon_mean[name] = carbon_global[name].collapsed('time', iris.analysis.MEAN)
    global_carbon_anomaly[name] = carbon_global[name]-global_carbon_mean[name]
    
    heat[name], heat_sumz[name], heat_global[name] = calc_ohc(temp[name], rhodzt[name], area[name])
    global_heat_mean[name] = heat_global[name].collapsed('time', iris.analysis.MEAN)
    global_heat_anomaly[name] = heat_global[name]-global_heat_mean[name]




names = {'carbon', 'heat'}
a1, a2, a3, a4, a5 = calc_zonal_regions(area['aredi_800'])

SO = {}
SMid = {}
NMid = {}
arctic = {}
tropics = {}
for name in names:
    if name == 'carbon':
        SO[name], SMid[name], NMid[name], arctic[name], tropics[name] =  calc_zonal_regions(carbon_sumz['aredi_800'])
    if name == 'heat':
        SO[name], SMid[name], NMid[name], arctic[name], tropics[name] =  calc_zonal_regions(heat_sumz['aredi_800'])

for name in names:
    # Sum in each region:
    SO[name] = SO[name]*a1
    SMid[name] = SMid[name]*a2
    NMid[name] = NMid[name]*a3
    arctic[name] = arctic[name]*a4
    tropics[name] = tropics[name]*a5

    SO[name] = SO[name].collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    SMid[name] = SMid[name].collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    NMid[name] = NMid[name].collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    tropics[name] = tropics[name].collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    arctic[name] = arctic[name].collapsed(['longitude', 'latitude'], iris.analysis.SUM)

    SO_mean = SO[name].collapsed('time', iris.analysis.MEAN)
    SMid_mean = SMid[name].collapsed('time', iris.analysis.MEAN)
    NMid_mean = NMid[name].collapsed('time', iris.analysis.MEAN)
    tropics_mean = tropics[name].collapsed('time', iris.analysis.MEAN)
    arctic_mean = arctic[name].collapsed('time', iris.analysis.MEAN)

    SO[name] = SO[name]-SO_mean
    SMid[name] = SMid[name] - SMid_mean
    NMid[name] = NMid[name] - NMid_mean
    tropics[name] = tropics[name] - tropics_mean
    arctic[name] = arctic[name] - arctic_mean

    # Detrend with polynomial
    if name == 'carbon':
        x = np.arange(0, 500, 1)

        p1 = np.poly1d(np.polyfit(x, SO[name].data, 4))
        SO[name] = SO[name] - p1(x)

        p1 = np.poly1d(np.polyfit(x, SMid[name].data, 4))
        SMid[name] = SMid[name] - p1(x)

        p1 = np.poly1d(np.polyfit(x, NMid[name].data, 4))
        NMid[name] = NMid[name] - p1(x)

        p1 = np.poly1d(np.polyfit(x, tropics[name].data, 4))
        tropics[name] = tropics[name] - p1(x)

        p1 = np.poly1d(np.polyfit(x, arctic[name].data, 4))
        arctic[name] = arctic[name] - p1(x)


SO_percent = {}
SMid_percent = {}
NMid_percent = {}
arctic_percent = {}
tropics_percent = {}

for name in names:
    if name == 'carbon':
        SO_percent[name] = SO[name]/global_carbon_anomaly['aredi_800']
        SMid_percent[name] = SMid[name]/global_carbon_anomaly['aredi_800']
        NMid_percent[name] = NMid[name]/global_carbon_anomaly['aredi_800']
        arctic_percent[name] = arctic[name]/global_carbon_anomaly['aredi_800']
        tropics_percent[name] = tropics[name]/global_carbon_anomaly['aredi_800']

    if name == 'heat':
        SO_percent[name] = SO[name]/global_heat_anomaly['aredi_800']
        SMid_percent[name] = SMid[name]/global_heat_anomaly['aredi_800']
        NMid_percent[name] = NMid[name]/global_heat_anomaly['aredi_800']
        arctic_percent[name] = arctic[name]/global_heat_anomaly['aredi_800']
        tropics_percent[name] = tropics[name]/global_heat_anomaly['aredi_800']


f, ((ax1), (ax3), (ax5), (ax7), (ax9)) = plt.subplots(5, 1, sharex='col', sharey='row')
ax1.plot(SO['carbon'].data/1e15)
ax2 = ax1.twinx()
ax2.plot(SO['heat'].data/1e22, color = 'b', ls = '--')
ax1.set_ylim([-4, 4])
ax2.set_ylim([-6, 6])
ax1.set_title('(a) Southern Ocean', fontsize = 10)

ax3.plot(SMid['carbon'].data/1e15, color = 'g')
ax4 = ax3.twinx()
ax4.plot(SMid['heat'].data/1e22, color = 'g', ls = '--')
ax3.set_ylim([-4, 4])
ax4.set_ylim([-6, 6])
ax3.set_title('(b) Southern Mid-Latitudes', fontsize = 10)

ax5.plot(NMid['carbon'].data/1e15, color = 'r')
ax6 = ax5.twinx()
ax6.plot(NMid['heat'].data/1e22, color = 'r', ls = '--')
ax5.set_ylim([-4, 4])
ax6.set_ylim([-6, 6])
ax5.set_title('(c) Northern Mid-Latitudes', fontsize = 10)

ax7.plot(arctic['carbon'].data/1e15, color = 'purple')
ax8 = ax7.twinx()
ax8.plot(arctic['heat'].data/1e22, color = 'purple', ls = '--')
ax7.set_ylim([-4, 4])
ax8.set_ylim([-6, 6])
ax7.set_title('(d) Arctic', fontsize = 10)

ax9.plot(tropics['carbon'].data/1e15, color = 'c')
ax10 = ax9.twinx()
ax10.plot(tropics['heat'].data/1e22, color = 'c', ls = '--')
ax9.set_ylim([-4, 4])
ax10.set_ylim([-6, 6])
ax9.set_title('(e) Tropics', fontsize = 10)



'''
fig = plt.figure()
plt.plot(SO_percent['heat'].data, color = 'b', ls = '--')
plt.plot(SMid_percent['heat'].data, color = 'b', ls = '--')
'''
plt.show()
