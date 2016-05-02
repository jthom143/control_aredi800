import iris
import iris.analysis.cartography
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np


from calc_OHC_OCC import calc_ohc
from calc_region import calc_six_regions, calc_zonal_regions


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


#SO, SPac, SAtl, NPac, NAtl, Arctic = calc_six_regions(temp['aredi_800'])
SO, Smid, Nmid, arctic, tropics  = calc_zonal_regions(temp['aredi_800'])

fig = plt.figure()
clevs = np.arange(0, 5, 1)

map = 'Paired'

lats = Nmid.coord('latitude').points
lons = Nmid.coord('longitude').points
data = Nmid[0, 0, :, :].data * 0 
ax = fig.add_axes([0.1,0.1,0.8,0.9], projection=ccrs.PlateCarree())
ax.coastlines()
ax.contourf(lons, lats, data,color = 'b')

lats = Smid.coord('latitude').points
lons = Smid.coord('longitude').points
data = Smid[0, 0, :, :].data * 0 + 1
ax.contourf(lons, lats, data,clevs, cmap = map)

lats = SO.coord('latitude').points
lons = SO.coord('longitude').points
data = SO[0, 0, :, :].data * 0 + 2
ax.contourf(lons, lats, data,clevs, cmap = map)

lats = arctic.coord('latitude').points
lons = arctic.coord('longitude').points
data = arctic[0, 0, :, :].data * 0 + 3
ax.contourf(lons, lats, data,clevs, cmap = map)

lats = tropics.coord('latitude').points
lons = tropics.coord('longitude').points
data = tropics[0,0,:,:].data *0 + 4
ax.contourf(lons, lats, data, clevs, cmap = map)

plt.show()
