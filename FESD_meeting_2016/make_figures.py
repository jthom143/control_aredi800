## Script to create figures for FESD 2016 presentation

import iris
import numpy as np
import iris.analysis.cartography
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys # access system routines

#sys.path.append('~/control_aredi800/python')
sys.path.append('/Users/jordanthomas/control_aredi800/python')

import colormaps as cmaps
from calc_region import calc_six_regions, calc_zonal_regions
from smooth import smooth

# Load data from '/home/jthom143/control_aredi800/FESD_meeting_2016/data/'
#PATH = '/home/jthom143/control_aredi800/FESD_meeting_2016/data/'
PATH = '/Users/jordanthomas/control_aredi800/FESD_meeting_2016/data/'

temp_ws       = iris.load_cube(PATH + 'temp_ws.nc')
mld_ws        = iris.load_cube(PATH + 'mld_ws.nc')
carbon        = iris.load_cube(PATH + 'carbon.nc')
heat          = iris.load_cube(PATH + 'heat.nc')
carbon_sumz   = iris.load_cube(PATH + 'carbon_sumz.nc')
heat_sumz     = iris.load_cube(PATH + 'heat_sumz.nc')
G             = iris.load_cube(PATH + 'G.nc')
carbon_global = iris.load_cube(PATH + 'carbon_global.nc')
carbon_SH     = iris.load_cube(PATH + 'carbon_SH.nc')
heat_global   = iris.load_cube(PATH + 'heat_global.nc')
heat_SH       = iris.load_cube(PATH + 'heat_SH.nc')
area = iris.load_cube(PATH+'area_t.nc')     # Tracer Cell Area

# Calculate carbon anomaly
global_carbon_mean = carbon_global.collapsed('time', iris.analysis.MEAN)
global_carbon_anomaly = carbon_global.data-global_carbon_mean.data

mean_SH = carbon_SH.collapsed('time', iris.analysis.MEAN)
carbon_SH_anomaly = carbon_SH - mean_SH

# Detrend Carbon with a polynomial
x = np.arange(0, 500, 1)
p = np.poly1d(np.polyfit(x, global_carbon_anomaly/1e15, 4))
p2 = np.poly1d(np.polyfit(x, carbon_SH_anomaly.data/1e15, 4))

global_carbon_detrend = global_carbon_anomaly/1e15 - p(x)
carbon_SH_detrend = carbon_SH_anomaly.data/1e15 - p2(x)

# Calculate heat anomaly
global_heat_mean = heat_global.collapsed('time', iris.analysis.MEAN)
global_heat_anomaly = heat_global.data-global_heat_mean.data

mean_SH = heat_SH.collapsed('time', iris.analysis.MEAN)
heat_SH_anomaly = heat_SH - mean_SH

# Detrend Heat with a polynomial
p = np.poly1d(np.polyfit(x, global_heat_anomaly/1e22, 4))
p2 = np.poly1d(np.polyfit(x, heat_SH_anomaly.data/1e22, 4))
x = np.arange(0, 500, 1)

global_heat_detrend = global_heat_anomaly/1e22 - p(x)
heat_SH_detrend = heat_SH_anomaly.data/1e22 - p2(x)

### FIGURE 1 #################################################################################################################################
fig = plt.figure(figsize=(7, 3))

clevs = np.arange(-1, 2.75, 0.25)
time = np.arange(0, 500, 1)

ax1 = fig.add_axes([0.2, 0.2, 0.6, 0.7])
im = ax1.contourf(time, temp_ws[:,:-6].coord('tcell pstar').points, (temp_ws[:,:-6].data).T , clevs,cmap = cmaps.viridis, extend = 'both')
ax1.plot(time, mld_ws.data, color = 'k')
plt.gca().invert_yaxis()
plt.ylabel('Depth [dbars]')
plt.xlabel('Time [years]')
plt.title('Subsurface Ocean Temperature', fontsize = 11)

cax = fig.add_axes([0.85,0.2,0.03,0.7])
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label(r'$^o$C')
plt.savefig('figures/temp_mld_depth.png', bbox_inches='tight')
##############################################################################################################################################

### FIGURE 2 #################################################################################################################################
fig = plt.figure(figsize=(7.5, 7.5))
ax1 = fig.add_axes([0.2, 0.6, 0.7, 0.3])
ax1.plot(global_carbon_detrend, color = 'k', lw = 1.5)
ax1.plot(carbon_SH_detrend, color = 'k', ls = '--', lw = 1.5)
plt.ylabel('Carbon Anomaly [PgC]')
plt.text(290, 3, 'Global mean = 37000 PgC', fontsize = 10)

ax2 = fig.add_axes([0.2, 0.2, 0.7, 0.3])
ax2.plot(global_heat_detrend, color = 'k', lw = 1.5)
ax2.plot(heat_SH_detrend, color = 'k', ls = '--', lw = 1.5)
plt.ylabel('Heat Anomaly [$10^{22}$J]')
plt.xlabel('Time [years]')
plt.text(290, 4, 'Global mean = 15 x 10$^{26}$ J', fontsize = 10)
plt.savefig('figures/occ_ohc_timeseries.png', bbox_inches='tight')
##############################################################################################################################################

### FIGURE 3 #################################################################################################################################
fig = plt.figure(figsize=(6.5, 5.5))
ax1 = fig.add_axes([0.15, 0.6, 0.75, 0.35])
ax1.plot(global_carbon_detrend, color = 'k', lw = 1.5)
ax1.plot(carbon_SH_detrend, color = 'k', ls = '--', lw = 1.5)
plt.ylabel('Carbon Anomaly [PgC]')
ax1.set_xlim([60, 160])

ax2 = fig.add_axes([0.15, 0.15, 0.75, 0.35])
ax2.plot(global_heat_detrend, color = 'k', lw = 1.5)
ax2.plot(heat_SH_detrend, color = 'k', ls = '--', lw = 1.5)
plt.ylabel('Heat Anomaly [$10^{22}$J]')
plt.xlabel('Time [years]')
ax2.set_xlim([60, 160])
plt.savefig('figures/occ_ohc_timeseries_period1a.png', bbox_inches='tight')

# Figure 3b
p1_start = 60
p1_end = 160
period1 = heat_SH_detrend[p1_start:p1_end]

max = np.argmax(period1)
min = np.argmin(period1)

min = p1_start + min
max = p1_start + max

fig = plt.figure(figsize=(6.5, 5.5))
ax1 = fig.add_axes([0.15, 0.6, 0.75, 0.35])
ax1.plot(global_carbon_detrend, color = 'k', lw = 1.5)
ax1.axvline(min, color = 'r', lw = 1.5)
ax1.axvline(max, color = 'b', lw = 1.5)
ax1.plot(carbon_SH_detrend, color = 'k', ls = '--', lw = 1.5)
plt.ylabel('Carbon Anomaly [PgC]')
ax1.set_xlim([60, 160])
plt.text(min+2, 2, 'Convective', color = 'r', fontsize = 10)
plt.text(max-23, -2, 'Non-convective', color = 'b', fontsize = 10)


ax2 = fig.add_axes([0.15, 0.15, 0.75, 0.35])
ax2.plot(global_heat_detrend, color = 'k', lw = 1.5)
ax2.plot(heat_SH_detrend, color = 'k', ls = '--', lw = 1.5)
ax2.axvline(min, color = 'r', lw = 1.5)
ax2.axvline(max, color = 'b', lw = 1.5)
plt.ylabel('Heat Anomaly [$10^{22}$J]')
plt.xlabel('Time [years]')
ax2.set_xlim([60, 160])
plt.savefig('figures/occ_ohc_timeseries_period1b.png', bbox_inches='tight')

##############################################################################################################################################

### FIGURE 4 #################################################################################################################################
data_carbon = carbon.collapsed('longitude', iris.analysis.MEAN)
data_carbon = data_carbon[min,:,:]-data_carbon[max,:,:]

data_heat   = heat.collapsed('longitude', iris.analysis.MEAN)
data_heat   = data_heat[min,:,:]-data_heat[max,:,:]

clevs_carbon = np.arange(-16, 18, 2)
clevs_heat = np.arange(-3, 3.5, 0.5)

fig = plt.figure(figsize=(12, 4))
ax1 = fig.add_axes([0.15, 0.25, 0.3, 0.6])
im = ax1.contourf(data_carbon.coord('latitude').points, data_carbon.coord('tcell pstar').points, data_carbon.data, clevs_carbon, cmap = 'RdBu_r', extend = 'both')
ax1.set_ylim([0, 2000])
ax1.invert_yaxis()
ax1.set_ylabel('Depth [dbars]')
ax1.set_xlabel('Latitude')
ax1.set_title('(a) Subsurface Carbon', fontsize = 11)

cax = fig.add_axes([0.47,0.25,0.02,0.6])
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label(r'Pg C')


ax2 = fig.add_axes([0.6, 0.25, 0.3, 0.6])
im2 = ax2.contourf(data_heat.coord('latitude').points, data_heat.coord('tcell pstar').points, data_heat.data/1e8, clevs_heat, cmap = 'RdBu_r', extend = 'both')
ax2.set_ylim([0, 2000])
ax2.invert_yaxis()
ax2.set_xlabel('Latitude')
ax2.set_title('(b) Subsurface Heat', fontsize = 11)

cax = fig.add_axes([0.92,0.25,0.02,0.6])
cb = fig.colorbar(im2, cax=cax, orientation='vertical')
cb.set_label(r'10$^8$J')
plt.savefig('figures/occ_ohc_subsurface.png', bbox_inches='tight')
##############################################################################################################################################

### FIGURE 5 #################################################################################################################################
names = {'carbon', 'heat'}
a1, a2, a3, a4, a5 = calc_zonal_regions(area)

SO = {}
SMid = {}
NMid = {}
arctic = {}
tropics = {}
SO_smooth = {}
SMid_smooth = {}
NMid_smooth = {}
arctic_smooth = {}
tropics_smooth = {}

for name in names:
    if name == 'carbon':
        SO[name], SMid[name], NMid[name], arctic[name], tropics[name] =  calc_zonal_regions(carbon_sumz)
    if name == 'heat':
        SO[name], SMid[name], NMid[name], arctic[name], tropics[name] =  calc_zonal_regions(heat_sumz)

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

    x = 21
    SO_smooth[name] = smooth(SO[name].data, window_len=x, window='hanning')
    SMid_smooth[name] = smooth(SMid[name].data, window_len=x, window='hanning')
    NMid_smooth[name] = smooth(NMid[name].data, window_len=x, window='hanning')
    tropics_smooth[name] = smooth(tropics[name].data, window_len=x, window='hanning')
    arctic_smooth[name] = smooth(arctic[name].data, window_len=x, window='hanning')

carbon_smooth = smooth(global_carbon_detrend, window_len=x, window='hanning')
heat_smooth = smooth(global_heat_detrend, window_len=x, window='hanning')


fig = plt.figure(figsize=(10, 4))
ax1 = fig.add_axes([0.1, 0.25, 0.37, 0.6])
ax1.plot(carbon_smooth, lw = 1.5, color = 'k')
ax1.plot(SO_smooth['carbon']/1e15, lw = 1.5, color = 'b')
ax1.plot(SMid_smooth['carbon']/1e15, lw = 1.5, color = 'g')
ax1.plot(NMid_smooth['carbon']/1e15, lw = 1.5, color = 'r')
ax1.plot(tropics_smooth['carbon']/1e15, lw = 1.5, color = 'purple')
ax1.plot(arctic_smooth['carbon']/1e15, lw = 1.5, color = 'c')

ax1.set_ylabel('Carbon Content [Pg C]')
ax1.set_xlabel('Time [years]')
ax1.set_title('(a) Carbon', fontsize = 11)
ax1.set_xlim([p1_start, p1_end])

plt.text(130, 3.6, 'Global', color = 'k', fontsize = 9)
plt.text(130, 3.2, 'Southern Ocean', color = 'b', fontsize = 9)
plt.text(130, 2.8, 'S. Mid-Lats', color = 'g', fontsize = 9)
plt.text(130, 2.4, 'N. Mid-Lats', color = 'r', fontsize = 9)
plt.text(130, 2.0, 'Arctic', color = 'c', fontsize = 9)
plt.text(130, 1.6, 'Tropics', color = 'purple', fontsize = 9)


ax2 = fig.add_axes([0.6, 0.25, 0.37, 0.6])
ax2.plot(heat_smooth, lw = 1.5, color = 'k')
ax2.plot(SO_smooth['heat']/1e22, lw = 1.5, color = 'b')
ax2.plot(SMid_smooth['heat']/1e22, lw = 1.5, color = 'g')
ax2.plot(NMid_smooth['heat']/1e22, lw = 1.5, color = 'r')
ax2.plot(tropics_smooth['heat']/1e22, lw = 1.5, color = 'purple')
ax2.plot(arctic_smooth['heat']/1e22, lw = 1.5, color = 'c')
ax2.set_xlabel('Time [years]')
ax2.set_ylabel('Heat Content [10$^{22}$J]')
ax2.set_title('(b) Heat', fontsize = 11)
ax2.set_xlim([p1_start, p1_end])

plt.savefig('figures/occ_ohc_regions.png', bbox_inches='tight')

##############################################################################################################################################

### FIGURE 6 #################################################################################################################################
f, ((ax1), (ax3), (ax5), (ax7), (ax9)) = plt.subplots(5, 1, sharex='col', sharey='row', figsize =(10, 10) )
ax1.plot(SO['carbon'].data/1e15, lw = 1.5)
ax2 = ax1.twinx()
ax2.plot(SO['heat'].data/1e22, color = 'b', ls = '--', lw = 1.5)
ax1.set_ylim([-4, 4])
ax2.set_ylim([-6, 6])
ax1.set_title('(a) Southern Ocean', fontsize = 11)

ax3.plot(SMid['carbon'].data/1e15, color = 'g', lw = 1.5)
ax4 = ax3.twinx()
ax4.plot(SMid['heat'].data/1e22, color = 'g', ls = '--', lw = 1.5)
ax3.set_ylim([-4, 4])
ax4.set_ylim([-6, 6])
ax3.set_title('(b) Southern Mid-Latitudes', fontsize = 11)

ax5.plot(NMid['carbon'].data/1e15, color = 'r', lw = 1.5)
ax6 = ax5.twinx()
ax6.plot(NMid['heat'].data/1e22, color = 'r', ls = '--', lw = 1.5)
ax5.set_ylim([-4, 4])
ax6.set_ylim([-6, 6])
ax5.set_title('(c) Northern Mid-Latitudes', fontsize = 11)
ax5.set_ylabel('Carbon Content [Pg C]')
ax6.set_ylabel('Heat Content [10 $^{22}$ J]')

ax7.plot(arctic['carbon'].data/1e15, color = 'c', lw = 1.5)
ax8 = ax7.twinx()
ax8.plot(arctic['heat'].data/1e22, color = 'c', ls = '--', lw = 1.5)
ax7.set_ylim([-4, 4])
ax8.set_ylim([-6, 6])
ax7.set_title('(d) Arctic', fontsize = 11)

ax9.plot(tropics['carbon'].data/1e15, color = 'purple', lw = 1.5)
ax10 = ax9.twinx()
ax10.plot(tropics['heat'].data/1e22, color = 'purple', ls = '--', lw = 1.5)
ax9.set_ylim([-4, 4])
ax10.set_ylim([-6, 6])
ax9.set_title('(e) Tropics', fontsize = 11)
plt.savefig('figures/compare_regions.png', bbox_inches='tight')


plt.show()
