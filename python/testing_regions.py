import iris
import iris.analysis.cartography
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np


import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_OHC_OCC import calc_ohc
from calc_region import calc_six_regions


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


SO, SPac, SAtl, NPac, NAtl, Arctic = calc_six_regions(temp['aredi_800'])

fig = plt.figure()
clevs = np.arange(1, 8, 1)

map = 'Paired'

lats = NAtl.coord('latitude').points
lons = NAtl.coord('longitude').points
data = NAtl[0, 0, :, :].data * 0 + 1
ax = fig.add_axes([0.1,0.1,0.8,0.9], projection=ccrs.PlateCarree())
ax.coastlines()
ax.contourf(lons, lats, data, clevs, cmap = map)

lats = SAtl.coord('latitude').points
lons = SAtl.coord('longitude').points
data = SAtl[0, 0, :, :].data * 0 + 2
ax.contourf(lons, lats, data,clevs, cmap = map)

lats = SPac.coord('latitude').points
lons = SPac.coord('longitude').points
data = SPac[0, 0, :, :].data * 0 + 3
ax.contourf(lons, lats, data,clevs, cmap = map)

lats = NPac.coord('latitude').points
lons = NPac.coord('longitude').points
data = NPac[0, 0, :, :].data * 0 + 4
ax.contourf(lons, lats, data,clevs, cmap = map)

lats = SO.coord('latitude').points
lons = SO.coord('longitude').points
data = SO[0, 0, :, :].data * 0 + 5
ax.contourf(lons, lats, data,clevs, cmap = map)

lats = Arctic.coord('latitude').points
lons = Arctic.coord('longitude').points
data = Arctic[0, 0, :, :].data * 0 + 6
ax.contourf(lons, lats, data,clevs, cmap = map)

plt.show()
