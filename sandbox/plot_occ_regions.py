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
from calc_region import calc_six_regions, calc_zonal_regions
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

names = ['aredi_800']

carbon = {}
carbon_sumz = {}
carbon_global = {}
global_carbon_anomaly = {}
global_carbon_mean = {}

for name in names:
    carbon[name], carbon_sumz[name], carbon_global[name] = calc_occ(dic[name], rhodzt[name], area[name])
    global_carbon_mean[name] = carbon_global[name].collapsed('time', iris.analysis.MEAN)
    global_carbon_anomaly[name] = carbon_global[name].data-global_carbon_mean[name].data

#SO, SPac, SAtl, NPac, NAtl, Arctic = calc_six_regions(carbon_sumz['aredi_800'])
#a1, a2, a3, a4, a5, a6 = calc_six_regions(area['aredi_800'])

SO, SMid, NMid, arctic, tropics =  calc_zonal_regions(carbon_sumz['aredi_800'])
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
plt.plot(SO.data/1e15, label = 'Southern Ocean')
plt.plot(SMid.data/1e15, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid.data/1e15, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics.data/1e15, label = 'Tropics')
plt.plot(arctic.data/1e15, label = 'Arctic')
plt.plot(global_carbon_anomaly['aredi_800'].data/1e15, color = 'k', lw = 2)
plt.legend()


# Detrend
SO_detrend = signal.detrend(SO.data)
SMid_detrend = signal.detrend(SMid.data)
NMid_detrend = signal.detrend(NMid.data)
tropics_detrend = signal.detrend(tropics.data)
arctic_detrend = signal.detrend(arctic.data) 
global_carbon_detrend = signal.detrend(global_carbon_anomaly['aredi_800'].data)

fig = plt.figure()
plt.plot(SO_detrend/1e15, label = 'Southern Ocean')
plt.plot(SMid_detrend/1e15, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid_detrend/1e15, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics_detrend/1e15, label = 'Tropics')
plt.plot(arctic_detrend/1e15, label = 'Arctic')
plt.plot(global_carbon_detrend/1e15, color = 'k', lw = 2)
plt.legend()
plt.xlabel('Time [years]')
plt.ylabel('Carbon content anomaly [Pg C]')

x = 21
SO_smooth = smooth(SO_detrend, window_len=x, window='hanning')
SMid_smooth = smooth(SMid_detrend, window_len=x, window='hanning')
NMid_smooth = smooth(NMid_detrend, window_len=x, window='hanning')
tropics_smooth = smooth(tropics_detrend, window_len=x, window='hanning')
arctic_smooth = smooth(arctic_detrend, window_len=x, window='hanning')
global_carbon_smooth = smooth(global_carbon_detrend, window_len = x, window='hanning')

fig = plt.figure(figsize=(12,8))
plt.plot(SO_smooth/1e15, label = 'Southern Ocean')
plt.plot(SMid_smooth/1e15, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid_smooth/1e15, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics_smooth/1e15, label = 'Tropics')
plt.plot(arctic_smooth/1e15, label = 'Arctic')
plt.plot(global_carbon_smooth/1e15, color = 'k', lw = 2)

plt.text(390,3.6, 'Southern Ocean', color = 'b')
plt.text(390,3.3, 'SH Mid-Lats', color = 'g')
plt.text(390,3.0, 'NH Mid-Lats', color = 'r')
plt.text(390,2.7, 'Tropics', color = 'c')
plt.text(390,2.4, 'Arctic', color = 'purple')
plt.text(390,2.1, 'Global', color = 'k')
plt.xlim([0, 500])
plt.xlabel('Time [years]')
plt.ylabel('Carbon content anomaly [Pg C]')
plt.savefig('/home/jthom143/control_aredi800/notes/figures/occ_regions.png')

fig = plt.figure()
sum = (SO_smooth + SMid_smooth +  NMid_smooth + tropics_smooth + arctic_smooth)/1e15
plt.plot(SO_smooth/1e15, label = 'Southern Ocean')
plt.plot(SMid_smooth/1e15, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid_smooth/1e15, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics_smooth/1e15, label = 'Tropics')
plt.plot(arctic_smooth/1e15, label = 'Arctic')

plt.plot(global_carbon_smooth/1e15, color = 'k', lw = 2)
plt.plot(sum, color = 'k', ls = '--', lw = 2)
plt.xlim([100, 200])
plt.xlabel('Time [years]')
plt.ylabel('Carbon content anomaly [Pg C]')


## Divide up into convective cycles
# Divide timeseries up into maxima/minima cycles
# Period 1: years 60-160
p1_start = 60
p1_end = 160

# Period 2: years 160-260
p2_start = 160
p2_end = 260

# Period 3: years 260-300
p3_start = 260
p3_end = 300

# Period 4: years 300-360
p4_start = 300
p4_end = 360

# Period 5: years 360-420
p5_start = 360
p5_end = 420

# Period 6: years 420-500
p6_start = 420


f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
ax1.plot(SO_detrend[p1_start:p1_end]/1e15, label = 'Southern Ocean')
ax1.plot(SMid_detrend[p1_start:p1_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax1.plot(NMid_detrend[p1_start:p1_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax1.plot(tropics_detrend[p1_start:p1_end]/1e15, label = 'Tropics')
ax1.plot(arctic_detrend[p1_start:p1_end]/1e15, label = 'Arctic')
ax1.plot(global_carbon_detrend[p1_start:p1_end]/1e15, color = 'k', lw = 2)
ax1.axhline(0, color = 'k', ls = '--')
ax1.set_title('Period 1 **', fontsize = 10)


ax2.plot(SO_detrend[p2_start:p2_end]/1e15, label = 'Southern Ocean')
ax2.plot(SMid_detrend[p2_start:p2_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax2.plot(NMid_detrend[p2_start:p2_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax2.plot(tropics_detrend[p2_start:p2_end]/1e15, label = 'Tropics')
ax2.plot(arctic_detrend[p2_start:p2_end]/1e15, label = 'Arctic')
ax2.plot(global_carbon_detrend[p2_start:p2_end]/1e15, color = 'k', lw = 2)
ax2.axhline(0, color = 'k', ls = '--')
ax2.set_title('Period 2 **', fontsize = 10)

ax3.plot(SO_detrend[p3_start:p3_end]/1e15, label = 'Southern Ocean')
ax3.plot(SMid_detrend[p3_start:p3_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax3.plot(NMid_detrend[p3_start:p3_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax3.plot(tropics_detrend[p3_start:p3_end]/1e15, label = 'Tropics')
ax3.plot(arctic_detrend[p3_start:p3_end]/1e15, label = 'Arctic')
ax3.plot(global_carbon_detrend[p3_start:p3_end]/1e15, color = 'k', lw = 2)
ax3.axhline(0, color = 'k', ls = '--')
ax3.set_title('Period 3', fontsize = 10)

ax4.plot(SO_detrend[p4_start:p4_end]/1e15, label = 'Southern Ocean')
ax4.plot(SMid_detrend[p4_start:p4_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax4.plot(NMid_detrend[p4_start:p4_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax4.plot(tropics_detrend[p4_start:p4_end]/1e15, label = 'Tropics')
ax4.plot(arctic_detrend[p4_start:p4_end]/1e15, label = 'Arctic')
ax4.plot(global_carbon_detrend[p4_start:p4_end]/1e15, color = 'k', lw = 2)
ax4.axhline(0, color = 'k', ls = '--')
ax4.set_title('Period 4', fontsize = 10)

ax5.plot(SO_detrend[p5_start:p5_end]/1e15, label = 'Southern Ocean')
ax5.plot(SMid_detrend[p5_start:p5_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax5.plot(NMid_detrend[p5_start:p5_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax5.plot(tropics_detrend[p5_start:p5_end]/1e15, label = 'Tropics')
ax5.plot(arctic_detrend[p5_start:p5_end]/1e15, label = 'Arctic')
ax5.plot(global_carbon_detrend[p5_start:p5_end]/1e15, color = 'k', lw = 2)
ax5.axhline(0, color = 'k', ls = '--')
ax5.set_title('Period 5 **', fontsize = 10)

ax6.plot(SO_detrend[p6_start:]/1e15, label = 'Southern Ocean')
ax6.plot(SMid_detrend[p6_start:]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax6.plot(NMid_detrend[p6_start:]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax6.plot(tropics_detrend[p6_start:]/1e15, label = 'Tropics')
ax6.plot(arctic_detrend[p6_start:]/1e15, label = 'Arctic')
ax6.plot(global_carbon_detrend[p6_start:]/1e15, color = 'k', lw = 2)
ax6.axhline(0, color = 'k', ls = '--')
ax6.set_title('Period 6', fontsize = 10)




### Detrended


f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
ax1.plot(SO_smooth[p1_start:p1_end]/1e15, label = 'Southern Ocean')
ax1.plot(SMid_smooth[p1_start:p1_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax1.plot(NMid_smooth[p1_start:p1_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax1.plot(tropics_smooth[p1_start:p1_end]/1e15, label = 'Tropics')
ax1.plot(arctic_smooth[p1_start:p1_end]/1e15, label = 'Arctic')
ax1.plot(global_carbon_smooth[p1_start:p1_end]/1e15, color = 'k', lw = 2)
ax1.axhline(0, color = 'k', ls = '--')
ax1.set_title('Period 1 **', fontsize = 10)


ax2.plot(SO_smooth[p2_start:p2_end]/1e15, label = 'Southern Ocean')
ax2.plot(SMid_smooth[p2_start:p2_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax2.plot(NMid_smooth[p2_start:p2_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax2.plot(tropics_smooth[p2_start:p2_end]/1e15, label = 'Tropics')
ax2.plot(arctic_smooth[p2_start:p2_end]/1e15, label = 'Arctic')
ax2.plot(global_carbon_smooth[p2_start:p2_end]/1e15, color = 'k', lw = 2)
ax2.axhline(0, color = 'k', ls = '--')
ax2.set_title('Period 2 **', fontsize = 10)

ax3.plot(SO_smooth[p3_start:p3_end]/1e15, label = 'Southern Ocean')
ax3.plot(SMid_smooth[p3_start:p3_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax3.plot(NMid_smooth[p3_start:p3_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax3.plot(tropics_smooth[p3_start:p3_end]/1e15, label = 'Tropics')
ax3.plot(arctic_smooth[p3_start:p3_end]/1e15, label = 'Arctic')
ax3.plot(global_carbon_smooth[p3_start:p3_end]/1e15, color = 'k', lw = 2)
ax3.axhline(0, color = 'k', ls = '--')
ax3.set_title('Period 3', fontsize = 10)

ax4.plot(SO_smooth[p4_start:p4_end]/1e15, label = 'Southern Ocean')
ax4.plot(SMid_smooth[p4_start:p4_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax4.plot(NMid_smooth[p4_start:p4_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax4.plot(tropics_smooth[p4_start:p4_end]/1e15, label = 'Tropics')
ax4.plot(arctic_smooth[p4_start:p4_end]/1e15, label = 'Arctic')
ax4.plot(global_carbon_smooth[p4_start:p4_end]/1e15, color = 'k', lw = 2)
ax4.axhline(0, color = 'k', ls = '--')
ax4.set_title('Period 4', fontsize = 10)

ax5.plot(SO_smooth[p5_start:p5_end]/1e15, label = 'Southern Ocean')
ax5.plot(SMid_smooth[p5_start:p5_end]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax5.plot(NMid_smooth[p5_start:p5_end]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax5.plot(tropics_smooth[p5_start:p5_end]/1e15, label = 'Tropics')
ax5.plot(arctic_smooth[p5_start:p5_end]/1e15, label = 'Arctic')
ax5.plot(global_carbon_smooth[p5_start:p5_end]/1e15, color = 'k', lw = 2)
ax5.axhline(0, color = 'k', ls = '--')
ax5.set_title('Period 5 **', fontsize = 10)

ax6.plot(SO_smooth[p6_start:]/1e15, label = 'Southern Ocean')
ax6.plot(SMid_smooth[p6_start:]/1e15, label = 'Southern Hemisphere Mid-Lats')
ax6.plot(NMid_smooth[p6_start:]/1e15, label = 'Northern Hemisphere Mid-Lats')
ax6.plot(tropics_smooth[p6_start:]/1e15, label = 'Tropics')
ax6.plot(arctic_smooth[p6_start:]/1e15, label = 'Arctic')
ax6.plot(global_carbon_smooth[p6_start:]/1e15, color = 'k', lw = 2)
ax6.axhline(0, color = 'k', ls = '--')
ax6.set_title('Period 6', fontsize = 10)
plt.savefig('/home/jthom143/control_aredi800/notes/figures/occ_regions_periods.png')



f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
ax1.plot(SO_smooth[p1_start:p1_end]/1e15, label = 'Southern Ocean')
#ax1.plot((SO_smooth[p1_start:p1_end]+ tropics_smooth[p1_start:p1_end])/1e15, color =  'b', ls = '--')
#ax1.plot((SO_smooth[p1_start:p1_end]+ SMid_smooth[p1_start:p1_end])/1e15, color =  'g', ls = '--')
ax1.plot((SO_smooth[p1_start:p1_end]+ tropics_smooth[p1_start:p1_end]+ SMid_smooth[p1_start:p1_end])/1e15, color =  'k', ls = '--')
ax1.plot(global_carbon_smooth[p1_start:p1_end]/1e15, color = 'k', lw = 2)
ax1.axhline(0, color = 'k', ls = '--')
ax1.set_title('Period 1', fontsize = 10)


ax2.plot(SO_smooth[p2_start:p2_end]/1e15, label = 'Southern Ocean')
#ax2.plot((SO_smooth[p2_start:p2_end]+ tropics_smooth[p2_start:p2_end])/1e15, color =  'b', ls = '--')
#ax2.plot((SO_smooth[p2_start:p2_end]+ SMid_smooth[p2_start:p2_end])/1e15, color =  'g', ls = '--')
ax2.plot((SO_smooth[p2_start:p2_end]+ tropics_smooth[p2_start:p2_end]+ SMid_smooth[p2_start:p2_end])/1e15,  color =  'k', ls = '--')
ax2.plot(global_carbon_smooth[p2_start:p2_end]/1e15, color = 'k', lw = 2)
ax2.axhline(0, color = 'k', ls = '--')
ax2.set_title('Period 2', fontsize = 10)

ax3.plot(SO_smooth[p3_start:p3_end]/1e15, label = 'Southern Ocean')
#ax3.plot((SO_smooth[p3_start:p3_end]+ tropics_smooth[p3_start:p3_end])/1e15, color =  'b', ls = '--')
#ax3.plot((SO_smooth[p3_start:p3_end]+ SMid_smooth[p3_start:p3_end])/1e15, color =  'g', ls = '--')
ax3.plot((SO_smooth[p3_start:p3_end]+ tropics_smooth[p3_start:p3_end]+ SMid_smooth[p3_start:p3_end])/1e15,  color =  'k', ls = '--')
ax3.plot(global_carbon_smooth[p3_start:p3_end]/1e15, color = 'k', lw = 2)
ax3.axhline(0, color = 'k', ls = '--')
ax3.set_title('Period 3', fontsize = 10)

ax4.plot(SO_smooth[p4_start:p4_end]/1e15, label = 'Southern Ocean')
#ax4.plot((SO_smooth[p4_start:p4_end]+ tropics_smooth[p4_start:p4_end])/1e15, color =  'b', ls = '--')
#ax4.plot((SO_smooth[p4_start:p4_end]+ SMid_smooth[p4_start:p4_end])/1e15, color =  'g', ls = '--')
ax4.plot((SO_smooth[p4_start:p4_end]+ tropics_smooth[p4_start:p4_end]+ SMid_smooth[p4_start:p4_end])/1e15,  color =  'k', ls = '--')
ax4.plot(global_carbon_smooth[p4_start:p4_end]/1e15, color = 'k', lw = 2)
ax4.axhline(0, color = 'k', ls = '--')
ax4.set_title('Period 4', fontsize = 10)

ax5.plot(SO_smooth[p5_start:p5_end]/1e15, label = 'Southern Ocean')
#ax5.plot((SO_smooth[p5_start:p5_end]+ tropics_smooth[p5_start:p5_end])/1e15, color =  'b', ls = '--')
#ax5.plot((SO_smooth[p5_start:p5_end]+ SMid_smooth[p5_start:p5_end])/1e15, color =  'g', ls = '--')
ax5.plot((SO_smooth[p5_start:p5_end]+ tropics_smooth[p5_start:p5_end]+ SMid_smooth[p5_start:p5_end])/1e15,  color =  'k', ls = '--')
ax5.plot(global_carbon_smooth[p5_start:p5_end]/1e15, color = 'k', lw = 2)
ax5.axhline(0, color = 'k', ls = '--')
ax5.set_title('Period 5', fontsize = 10)

ax6.plot(SO_smooth[p6_start:]/1e15, label = 'Southern Ocean')
#ax6.plot((SO_smooth[p6_start:]+ tropics_smooth[p6_start:])/1e15, color =  'b', ls = '--')
#ax6.plot((SO_smooth[p6_start:]+ SMid_smooth[p6_start:])/1e15, color =  'g', ls = '--')
ax6.plot((SO_smooth[p6_start:]+ tropics_smooth[p6_start:]+ SMid_smooth[p6_start:])/1e15,  color =  'k', ls = '--')
ax6.plot(global_carbon_smooth[p6_start:]/1e15, color = 'k', lw = 2)
ax6.axhline(0, color = 'k', ls = '--')
ax6.set_title('Period 6', fontsize = 10)






plt.show()

