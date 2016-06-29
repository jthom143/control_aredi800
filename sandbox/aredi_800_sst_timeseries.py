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
import iris.analysis.cartography

import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_region import calc_zonal_regions


## Load Data ###
'''
    Files loaded in this script were created in 'plot_aredi_climatologies.py'
    '''

# Load Temp Data
PATHS={
       'aredi_800':'~/control_aredi800/data/newCO2_control_800/'
        }


sst= {}
area = {}


print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    sst[name] = iris.load_cube(PATH+'sst.nc')
    area[name] = iris.load_cube(PATH+'area_t.nc')
    
    if name == 'aredi_800':
        sst[name] = sst[name][-500:]
    if name == 'aredi_2400':
        sst[name] = sst[name][:500]
    if name == 'low_gm':
        sst[name] = sst[name][:500]

SO, SMid, NMid, arctic, tropics =  calc_zonal_regions(sst['aredi_800'])
a1, a2, a3, a4, a5 = calc_zonal_regions(area['aredi_800'])

# MEAN in each region:
SO.coord('latitude').guess_bounds()
SO.coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(SO)
SO = SO.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

SMid.coord('latitude').guess_bounds()
SMid.coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(SMid)
SMid = SMid.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

NMid.coord('latitude').guess_bounds()
NMid.coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(NMid)
NMid = NMid.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

tropics.coord('latitude').guess_bounds()
tropics.coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(tropics)
tropics = tropics.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

arctic.coord('latitude').guess_bounds()
arctic.coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(arctic)
arctic = arctic.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

sst['aredi_800'].coord('latitude').guess_bounds()
sst['aredi_800'].coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(sst['aredi_800'])
global_sst = sst['aredi_800'].collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights = grid_areas)

fig = plt.figure()
plt.plot(global_sst.data, color = 'k', label = 'Global')
plt.plot(SO.data, label = 'Southern Ocean')
plt.plot(SMid.data, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid.data, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics.data, label = 'Tropics')
plt.plot(arctic.data, label = 'Arctic')


# Calculate anomalies
SO_mean = SO.collapsed('time', iris.analysis.MEAN)
SMid_mean = SMid.collapsed('time', iris.analysis.MEAN)
NMid_mean = NMid.collapsed('time', iris.analysis.MEAN)
tropics_mean = tropics.collapsed('time', iris.analysis.MEAN)
arctic_mean = arctic.collapsed('time', iris.analysis.MEAN)
global_mean = global_sst.collapsed('time', iris.analysis.MEAN)

SO = SO-SO_mean
SMid = SMid - SMid_mean
NMid = NMid - NMid_mean
tropics = tropics - tropics_mean
arctic = arctic - arctic_mean
global_sst = global_sst - global_mean

fig = plt.figure()
plt.plot(SO.data, label = 'Southern Ocean')
plt.plot(SMid.data, label = 'Southern Hemisphere Mid-Lats')
plt.plot(NMid.data, label = 'Northern Hemisphere Mid-Lats')
plt.plot(tropics.data, label = 'Tropics')
plt.plot(arctic.data, label = 'Arctic')
plt.plot(global_sst.data, color = 'k', label = 'Global', lw = 2)
plt.legend()


f, ((ax1), (ax2), (ax3), (ax4), (ax5)) = plt.subplots(5, 1, sharex='col', sharey='row', figsize=(12,8))
ax1.plot(global_sst.data, color = 'k', label = 'Global', lw = 2)
ax1.plot(SO.data, color = 'b')
ax1.set_title('(a) Southern Ocean', fontsize = 10)
    
ax2.plot(global_sst.data, color = 'k', label = 'Global', lw = 2)
ax2.plot(SMid.data, color = 'g')
ax2.set_title('(b) S. Mid-Lats', fontsize = 10)

ax3.plot(global_sst.data, color = 'k', label = 'Global', lw = 2)
ax3.plot(NMid.data, color = 'r')
ax3.set_title('(c) N. Mid-Lats', fontsize = 10)
ax3.set_ylabel('SST [$^o$C]', fontsize = 10)

ax4.plot(global_sst.data, color = 'k', label = 'Global', lw = 2)
ax4.plot(tropics.data, color = 'c')
ax4.set_title('(d) Tropics', fontsize = 10)

ax5.plot(global_sst.data, color = 'k', label = 'Global', lw = 2)
ax5.plot(arctic.data, color = 'purple')
ax5.set_title('(e) Arctic', fontsize = 10)
plt.savefig('notes/figures/regional_ssts.png')



## Plot SST Contour for convective-non-convective
p1_min = 136
p1_max = 99

p2_min = 337
p2_max = 204

p5_min = 405
p5_max = 380

sst_p1_contour = sst['aredi_800'][p1_min,:,:] - sst['aredi_800'][p1_max,:,:]
sst_p2_contour = sst['aredi_800'][p2_min,:,:] - sst['aredi_800'][p2_max,:,:]
sst_p5_contour = sst['aredi_800'][p5_min,:,:] - sst['aredi_800'][p5_max,:,:]

lats = sst['aredi_800'].coord('latitude').points
lons = sst['aredi_800'].coord('longitude').points


fig = plt.figure()
clevs = np.arange(-3, 3.2, 0.2)
ax1 = plt.axes(projection=ccrs.PlateCarree())
im = ax1.contourf(lons, lats, sst_p1_contour.data, clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('$^{o}$C')
ax1.coastlines()
plt.title('(a) Period 1', fontsize = 11)
plt.savefig('notes/figures/sst_contour_1.png', bbox_inches='tight')



fig = plt.figure()
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.contourf(lons, lats, sst_p2_contour.data,clevs, cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('$^{o}$C')
ax2.coastlines()
plt.title('(b) Period 2', fontsize = 11)
plt.savefig('notes/figures/sst_contour_2.png', bbox_inches='tight')


fig = plt.figure()
ax3 = plt.axes(projection=ccrs.PlateCarree())
ax3.contourf(lons, lats, sst_p5_contour.data, clevs,cmap = 'RdBu_r', extend = 'both', transform=ccrs.PlateCarree())
cb = plt.colorbar(im, orientation='horizontal')
cb.set_label('$^{o}$C')
ax3.coastlines()
plt.title('(c) Period 3', fontsize = 11)
plt.savefig('notes/figures/sst_contour_3.png', bbox_inches='tight')

iris.save(sst_p1_contour, 'data/sst_p1_contour.nc')
iris.save(sst_p2_contour, 'data/sst_p2_contour.nc')
iris.save(sst_p5_contour, 'data/sst_p5_contour.nc')

iris.save(sst['aredi_800'], 'data/sst_cube.nc')


############################
## Calculate sigma*T^4 to use in GHG analysis
sigma = 5.67e-8
sst_K = sst['aredi_800'] + 273.15
stephan_boltzmann = sigma*(sst_K**4)
iris.save(stephan_boltzmann, 'data/stephan_boltzmann.nc')

plt.show()

