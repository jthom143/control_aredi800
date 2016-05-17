
"""
    Script to recreate figure 1 from Bernadello et al, 2004
    """

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
import glob
import time
import numpy as np
import cartopy.crs as ccrs
import matplotlib
import cartopy.feature
import matplotlib.path as mpath
import matplotlib.patches as mpatches

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

    
    if name == 'aredi_800':
        temp[name] = temp[name][-500:]
    if name == 'aredi_2400':
        temp[name] = temp[name][:500]
        rhodzt[name] = rhodzt[name][:500]
        mld[name] = mld[name][:500]
    if name == 'low_gm':
        temp[name] = temp[name][:500]
        rhodzt[name] = rhodzt[name][:500]
        mld[name] = mld[name][:500]


## Calculate OHC ###
names = {'aredi_400', 'aredi_800', 'aredi_2400', 'low_gm'}

temp_ws = {}
mld_ws = {}

for name in names:
    temp_ws[name]= calc_weddell_sea(temp[name])
    mld_ws[name] = calc_weddell_sea(mld[name])
    
    mld_ws[name].coord('latitude').guess_bounds()
    mld_ws[name].coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mld_ws[name])
    mld_ws[name] = mld_ws[name].collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)

    temp_ws[name].coord('latitude').guess_bounds()
    temp_ws[name].coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(temp_ws[name])
    temp_ws[name] = temp_ws[name].collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)



clevs = np.arange(-1, 2.75, 0.25)
time = np.arange(0, 500, 1)


### Create Figure
fig = plt.figure(figsize=(6, 8))
ax2 = fig.add_axes([0.1, 0.79, 0.7, 0.15])
ax2.contourf(time, temp_ws['aredi_400'][:,:-6].coord('tcell pstar').points, (temp_ws['aredi_400'][:,:-6].data).T, clevs,cmap = cmaps.viridis, extend = 'both')
ax2.plot(time, mld_ws['aredi_400'].data, color = 'k')
plt.gca().invert_yaxis()
#plt.xticks([2000,2100, 2200, 2300, 2400, 2500], [' ',' ', ' ', ' ', ' ', ' '])
plt.title(r'(a) A$_{\rm redi}$ = 400 m$^2$s$^{-1}$', fontsize = 11)


ax3 = fig.add_axes([0.1, 0.56, 0.7, 0.15])
ax3.contourf(time, temp_ws['aredi_800'][:,:-6].coord('tcell pstar').points, (temp_ws['aredi_800'][:,:-6].data).T , clevs,cmap = cmaps.viridis, extend = 'both')
ax3.plot(time, mld_ws['aredi_800'].data, color = 'k')
plt.gca().invert_yaxis()
plt.ylabel('Depth (dbars)')
#plt.xticks([2000,2100, 2200, 2300, 2400, 2500], [' ',' ', ' ', ' ', ' ', ' '])
plt.title(r'(b) A$_{\rm redi}$ = 800 m$^2$s$^{-1}$', fontsize = 11)

ax4 = fig.add_axes([0.1, 0.33, 0.7, 0.15])
im2 = ax4.contourf(time, temp_ws['aredi_2400'][:,:-6].coord('tcell pstar').points, (temp_ws['aredi_2400'][:,:-6].data).T , clevs,cmap =  cmaps.viridis, extend = 'both')
ax4.plot(time, mld_ws['aredi_2400'].data, color = 'k')
plt.gca().invert_yaxis()
plt.title(r'(c) A$_{\rm redi}$ = 2400 m$^2$s$^{-1}$', fontsize = 11)

ax4 = fig.add_axes([0.1, 0.1, 0.7, 0.15])
im2 = ax4.contourf(time, temp_ws['low_gm'][:,:-6].coord('tcell pstar').points, (temp_ws['low_gm'][:,:-6].data).T , clevs,cmap =  cmaps.viridis, extend = 'both')
ax4.plot(time, mld_ws['low_gm'].data, color = 'k')
plt.gca().invert_yaxis()
plt.title(r'(d) GM$_{\rm min}$ = 600 m$^2$s$^{-1}$', fontsize = 11)
plt.xlabel('Time (years)')

cax = fig.add_axes([0.85,0.3,0.03,0.4])
cb = fig.colorbar(im2, cax=cax, orientation='vertical')
cb.set_label('$^o$C')
plt.savefig('notes/figures/heat_depth.png', bbox_inches='tight')


plt.show()
