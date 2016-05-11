"""
    Script to plot occ scatter plots for different aredi simulations
"""

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import cartopy.crs as ccrs
import matplotlib
from netCDF4 import Dataset
from scipy import signal

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_occ
from calc_convective_periods import convect

## Load Data ###
'''
   Files loaded in this script were created in 'plot_aredi_climatologies.py'
'''

# Load Temp Data
PATHS={'aredi_400':'~/control_aredi800/data/newCO2_control_400_ag/',
       'aredi_800':'~/control_aredi800/data/newCO2_control_800/',
       'aredi_2400':'~/control_aredi800/data/newCO2_control_2400_ag/'
       'low_gm': '~/control_aredi800/data/gm_min_600/'
              }

dic = {}
rhodzt = {}
area = {}
mld = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
	print name
        dic[name] = iris.load_cube(PATH+'dic.nc')
        rhodzt[name] = iris.load_cube(PATH+'rho_dzt.nc')
        area[name] = iris.load_cube(PATH+'area_t.nc')
	mld[name] = iris.load_cube(PATH+'mld.nc')

	if name == 'aredi_800':
	    dic[name] = dic[name][-500:]
	    
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
	global_carbon_detrend[name] = signal.detrend(global_carbon_anomaly[name])

	# Southern Hemisphere:
	constraint = iris.Constraint(latitude=lambda y: -90 < y < -50)
	carbon_SH_tmp = carbon_sumz[name].extract(constraint)
	area_SH = area[name].extract(constraint)
	b = carbon_SH_tmp*area_SH
	carbon_SH[name] = b.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
	mean_SH = carbon_SH[name].collapsed('time', iris.analysis.MEAN)
	carbon_SH_anomaly[name] = carbon_SH[name] - mean_SH
	carbon_SH_detrend[name] = signal.detrend(carbon_SH_anomaly[name].data)

## Calculate convective years ##
year_convect = {}

for name in names:
	year_convect[name], mld_convect, area_convect = convect(mld[name], area[name])

# Figure
fig1 = plt.figure(figsize = (5,10))
ax1 = fig1.add_axes([0.15, 0.75, 0.75, 0.15])
ax1.scatter((global_carbon_detrend['aredi_400'])/1e15, carbon_SH_detrend['aredi_400']/1e15, marker ='.', color = 'k') 
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11) 

ax2 = fig1.add_axes([0.15, 0.5, 0.75, 0.15])
ax2.scatter((global_carbon_detrend['aredi_800'])/1e15, carbon_SH_detrend['aredi_800']/1e15, marker ='.', color = 'k')
ax2.set_ylabel('SH Carbon Content Anomaly [Pg C]')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)

ax3 = fig1.add_axes([0.15, 0.25, 0.75, 0.15])
ax3.scatter((global_carbon_detrend['aredi_2400'])/1e15, carbon_SH_detrend['aredi_2400']/1e15, marker ='.', color = 'k')
ax3.set_xlabel('Global Carbon Content Anomaly [Pg C]')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
plt.savefig('notes/aredi_occ_scatter.png', bbox_inches='tight')

ax4 = fig1.add_axes([0.15, 0.1, 0.75, 0.15])
ax4.scatter((global_carbon_detrend['aredi_2400'])/1e15, carbon_SH_detrend['aredi_2400']/1e15, marker ='.', color = 'k')
ax4.set_xlabel('Global Carbon Content Anomaly [Pg C]')
plt.title('(c) GM$_{min}$ = 600 m$^2$s$^{-1}$', fontsize = 11)
plt.savefig('notes/aredi_occ_scatter.png', bbox_inches='tight')


# Figure - highlight convective times
a_400 = ma.masked_where(ma.getmask(year_convect['aredi_400']), global_carbon_detrend['aredi_400'])
b_400 = ma.masked_where(ma.getmask(year_convect['aredi_400']), carbon_SH_detrend['aredi_400'])

a_800 = ma.masked_where(ma.getmask(year_convect['aredi_800']), global_carbon_detrend['aredi_800'])
b_800 = ma.masked_where(ma.getmask(year_convect['aredi_800']), carbon_SH_detrend['aredi_800'])

a_2400 = ma.masked_where(ma.getmask(year_convect['aredi_2400']), global_carbon_detrend['aredi_2400'])
b_2400 = ma.masked_where(ma.getmask(year_convect['aredi_2400']), carbon_SH_detrend['aredi_2400'])

fig1 = plt.figure(figsize = (5,10))
ax1 = fig1.add_axes([0.15, 0.75, 0.75, 0.15])
ax1.scatter((global_carbon_detrend['aredi_400'])/1e15, carbon_SH_detrend['aredi_400']/1e15, marker ='.', color = 'k') 
ax1.scatter(a_400/1e15, b_400/1e15, marker = '.', color = '#368BC1')
plt.title('(a) A$_{redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11) 

ax2 = fig1.add_axes([0.15, 0.5, 0.75, 0.15])
ax2.scatter((global_carbon_detrend['aredi_800'])/1e15, carbon_SH_detrend['aredi_800']/1e15, marker ='.', color = 'k')
ax2.scatter(a_800/1e15, b_800/1e15, marker = '.', color = '#368BC1')
ax2.set_ylabel('SH Carbon Content Anomaly [Pg C]')
plt.title('(b) A$_{redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)

ax3 = fig1.add_axes([0.15, 0.25, 0.75, 0.15])
ax3.scatter((global_carbon_detrend['aredi_2400'])/1e15, carbon_SH_detrend['aredi_2400']/1e15, marker ='.', color = 'k')
ax3.scatter(a_2400/1e15, b_2400/1e15, marker = '.', color = '#368BC1')
ax3.set_xlabel('Global Carbon Content Anomaly [Pg C]')
plt.title('(c) A$_{redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)
plt.savefig('notes/aredi_occ_scatter_wconv.png', bbox_inches='tight')





plt.show()
