"""
    Script to plot the timeseries of the aredi 800 simulation
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
import sys # access system routines
sys.path.append('/home/jthom143/python_functions')
import colormaps as cmaps


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



# Divide timeseries up into maxima/minima cycles
# Period 1: years 60-160
p1_start = 60
p1_end = 160
period1 = global_carbon_detrend['aredi_800'][p1_start:p1_end]
period1_sh = carbon_SH_detrend['aredi_800'][p1_start:p1_end]

max = np.argmax(period1)
max_sh = np.argmax(period1_sh)

min = np.argmin(period1)
min_sh = np.argmin(period1_sh)

# Period 2: years 160-260
p2_start = 160
p2_end = 260
period2 = global_carbon_detrend['aredi_800'][p2_start:p2_end]
period2_sh = carbon_SH_detrend['aredi_800'][p2_start:p2_end]

max2 = np.argmax(period2)
max2_sh = np.argmax(period2_sh)

min2 = np.argmin(period2)
min2_sh = np.argmin(period2_sh)

# Period 3: years 260-300
p3_start = 260
p3_end = 300
period3 = global_carbon_detrend['aredi_800'][p3_start:p3_end]
period3_sh = carbon_SH_detrend['aredi_800'][p3_start:p3_end]

max3 = np.argmax(period3)
max3_sh = np.argmax(period3_sh)

min3 = np.argmin(period3)
min3_sh = np.argmin(period3_sh)

# Period 4: years 300-360
p4_start = 300
p4_end = 360
period4 = global_carbon_detrend['aredi_800'][p4_start:p4_end]
period4_sh = carbon_SH_detrend['aredi_800'][p4_start:p4_end]

max4 = np.argmax(period4)
max4_sh = np.argmax(period4_sh)

min4 = np.argmin(period4)
min4_sh = np.argmin(period4_sh)

# Period 5: years 360-420
p5_start = 360
p5_end = 420
period5 = global_carbon_detrend['aredi_800'][p5_start:p5_end]
period5_sh = carbon_SH_detrend['aredi_800'][p5_start:p5_end]

max5 = np.argmax(period5)
max5_sh = np.argmax(period5_sh)

min5 = np.argmin(period5)
min5_sh = np.argmin(period5_sh)

# Period 6: years 420-500
p6_start = 420
period6 = global_carbon_detrend['aredi_800'][p6_start:]
period6_sh = carbon_SH_detrend['aredi_800'][p6_start:]

max6 = np.argmax(period6)
max6_sh = np.argmax(period6_sh)

min6 = np.argmin(period6)
min6_sh = np.argmin(period6_sh)

## For each Convective cycle, plot the ocean carbon depth vs latitude at each maximum:
f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p1_start+min_sh,:,:]-datas[p1_start+min,:,:]
clevs = np.arange(-16, 18, 2)

ax1.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')
ax1.invert_yaxis()


datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p2_start+min2_sh,:,:]-datas[p2_start+min2,:,:]
ax2.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')


datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p3_start+min3_sh,:,:]-datas[p3_start+min3,:,:]
ax3.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')
ax3.invert_yaxis()


datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p4_start+min4_sh,:,:]-datas[p4_start+min4,:,:]
ax4.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')


datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p5_start+min5_sh,:,:]-datas[p5_start+min5,:,:]
ax5.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')
ax5.invert_yaxis()


datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p6_start+min6_sh,:,:]-datas[p6_start+min6,:,:]
ax6.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')
plt.savefig('notes/figures/aredi800_occ_max_contour.png')





f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p1_start+min,:,:]-datas[p1_start+max,:,:]
clevs = np.arange(-16, 18, 2)

ax1.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')
ax1.invert_yaxis()

datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p2_start+min2,:,:]-datas[p2_start+max2,:,:]
ax2.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')


datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p3_start+min3,:,:]-datas[p3_start+max3,:,:]
ax3.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')
ax3.invert_yaxis()

datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p4_start+min4,:,:]-datas[p4_start+max4,:,:]
ax4.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')


datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p5_start+min5,:,:]-datas[p5_start+max5,:,:]
ax5.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')
ax5.invert_yaxis()

datas = carbon['aredi_800'].collapsed('longitude', iris.analysis.MEAN)
datas = datas[p6_start+min6_sh,:,:]-datas[p6_start+min6,:,:]
ax6.contourf(datas.coord('latitude').points, datas.coord('tcell pstar').points, datas.data, clevs, cmap = 'RdBu_r', extend = 'both')
plt.savefig('notes/figures/aredi800_occ_min-max_contour.png')





plt.show()




