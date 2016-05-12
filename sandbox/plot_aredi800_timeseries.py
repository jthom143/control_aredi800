  
"""
    Script to plot occ variance
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
    if name == 'aredi_2400':
        dic[name] = dic[name][:500]
        rhodzt[name] = rhodzt[name][:500]
    if name == 'low_gm':
        dic[name] = dic[name][:500]
        rhodzt[name] = rhodzt[name][:500]

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
    
    # Southern Hemisphere:
    constraint = iris.Constraint(latitude=lambda y: -90 < y < -50)
    carbon_SH_tmp = carbon_sumz[name].extract(constraint)
    area_SH = area[name].extract(constraint)
    b = carbon_SH_tmp*area_SH
    carbon_SH[name] = b.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    mean_SH = carbon_SH[name].collapsed('time', iris.analysis.MEAN)
    carbon_SH_anomaly[name] = carbon_SH[name] - mean_SH
    
    # Detrend with polynomial
    x = np.arange(0, 500, 1)
    
    p1 = np.poly1d(np.polyfit(x, global_carbon_anomaly[name]/1e15, 4))
    global_carbon_detrend[name] = global_carbon_anomaly[name]/1e15 - p1(x)
    
    p2 = np.poly1d(np.polyfit(x, carbon_SH_anomaly[name].data/1e15, 4))
    carbon_SH_detrend[name] = carbon_SH_anomaly[name].data/1e15 - p2(x)

fig1 = plt.figure(figsize = (8, 4))
plt.plot(global_carbon_detrend['aredi_800'], color = 'k', lw = 2)
plt.plot(carbon_SH_detrend['aredi_800'], color = 'k', lw = 2, ls = '--')
plt.savefig('notes/figures/aredi800_occ_timeseries.png')


# Divide timeseries up into maxima/minima cycles
# Period 1: years 60-160
period1 = global_carbon_detrend['aredi_800'][60:161]
period1_sh = carbon_SH_detrend['aredi_800'][60:161]

max = np.argmax(period1)
max_sh = np.argmax(period1_sh)

# Period 2: years 160-260
period2 = global_carbon_detrend['aredi_800'][160:261]
period2_sh = carbon_SH_detrend['aredi_800'][160:261]

max2 = np.argmax(period2)
max2_sh = np.argmax(period2_sh)

# Period 3: years 260-300
period3 = global_carbon_detrend['aredi_800'][260:301]
period3_sh = carbon_SH_detrend['aredi_800'][260:301]

max3 = np.argmax(period3)
max3_sh = np.argmax(period3_sh)

# Period 4: years 300-360
period4 = global_carbon_detrend['aredi_800'][300:361]
period4_sh = carbon_SH_detrend['aredi_800'][300:361]

max4 = np.argmax(period4)
max4_sh = np.argmax(period4_sh)

# Period 5: years 360-420
period5 = global_carbon_detrend['aredi_800'][360:421]
period5_sh = carbon_SH_detrend['aredi_800'][360:421]

max5 = np.argmax(period5)
max5_sh = np.argmax(period5_sh)

# Period 6: years 420-500
period6 = global_carbon_detrend['aredi_800'][420:]
period6_sh = carbon_SH_detrend['aredi_800'][420:]

max6 = np.argmax(period6)
max6_sh = np.argmax(period6_sh)

#fig3 = plt.figure(figsize = (8, 10))
f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
ax1.plot(period1, color = 'k', lw = 1.5)
ax1.plot(period1_sh, color = 'k', lw = 1.5, ls = '--')
ax1.axvline(max, color = 'b', lw = 1.5)
ax1.axvline(max_sh, color = 'b', lw = 1.5, ls = '--')

ax2.plot(period2, color = 'k', lw = 1.5)
ax2.plot(period2_sh, color = 'k', lw = 1.5, ls = '--')
ax2.axvline(max2, color = 'b', lw = 1.5)
ax2.axvline(max2_sh, color = 'b', lw = 1.5, ls = '--')

ax3.plot(period3, color = 'k', lw = 1.5)
ax3.plot(period3_sh, color = 'k', lw = 1.5, ls = '--')
ax3.axvline(max3, color = 'b', lw = 1.5)
ax3.axvline(max3_sh, color = 'b', lw = 1.5, ls = '--')

ax4.plot(period4, color = 'k', lw = 1.5)
ax4.plot(period4_sh, color = 'k', lw = 1.5, ls = '--')
ax4.axvline(max4, color = 'b', lw = 1.5)
ax4.axvline(max4_sh, color = 'b', lw = 1.5, ls = '--')

ax5.plot(period5, color = 'k', lw = 1.5)
ax5.plot(period5_sh, color = 'k', lw = 1.5, ls = '--')
ax5.axvline(max5, color = 'b', lw = 1.5)
ax5.axvline(max5_sh, color = 'b', lw = 1.5, ls = '--')

ax6.plot(period6, color = 'k', lw = 1.5)
ax6.plot(period6_sh, color = 'k', lw = 1.5, ls = '--')
ax6.axvline(max6, color = 'b', lw = 1.5)
ax6.axvline(max6_sh, color = 'b', lw = 1.5, ls = '--')
plt.savefig('notes/figures/aredi800_occ_timelag.png')



print max_sh-max
print max2_sh-max2
print max3_sh-max3
print max4_sh-max4
print max5_sh-max5
print max6_sh-max6







plt.show()




