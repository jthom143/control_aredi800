"""
Code to define convective time-periods in GFDL ESM2Mc aredi=800 control run
"""

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
import glob
import time
import cartopy.crs as ccrs
import cartopy.feature
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import numpy.ma as ma
import matplotlib.collections as collections

# Import my functions
import sys # access system routines
sys.path.append('/home/jthom143/python_functions')
from calc_region import calc_weddell_sea 


def convect(mld, area, mld_avg):
    import numpy as np
    import iris
    import numpy.ma as ma

    t = len(mld.coord('time').points)

    # Create area array with time component
    A = area.data
    area_ts = np.tile(A[np.newaxis, :,:], (t,1,1))

    # Extract when MLD is greater than 1000m
    mld_600 = ma.masked_where(mld.data< 1000, mld.data)
    area_600 = ma.masked_where(mld.data< 1000, area_ts)

    # Convert to Cube
    mld_600 = iris.cube.Cube(mld_600, long_name = 'Convective mld', units = 'm', dim_coords_and_dims = [(mld.coord('time'), 0), (mld.coord('latitude'), 1), (mld.coord('longitude'), 2)])
    area_convection = iris.cube.Cube(area_600, long_name = 'Convective area', units = 'm^2', dim_coords_and_dims = [(mld.coord('time'), 0), (area.coord('latitude'), 1), (area.coord('longitude'), 2)])
    area_convection = area_convection.collapsed(['latitude', 'longitude'], iris.analysis.SUM)
    area_convection.convert_units('km^2')

    # Extract when convection area > 100,000 km^2

    year_convect = ma.masked_where(area_convection.data<100000, mld.coord('time').points)
    mld_convect = ma.masked_where(area_convection.data<100000, mld_avg.data)
    area_convect = ma.masked_where(area_convection.data<100000, area_convection.data)

    return year_convect, mld_convect, area_convect
                                                                            



# Load Data
PATH = '~/control_aredi800/data/'

mld = iris.load_cube(PATH+'mld.nc')
area = iris.load_cube(PATH+'area_t.nc')

## Isolate Weddell Sea
mld_ws = calc_weddell_sea(mld)
area_ws = calc_weddell_sea(area)

mld_ws.coord('latitude').guess_bounds()
mld_ws.coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(mld_ws)
total_mld_ws = mld_ws.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)

# Calculate convection metrics
year_convect, mld_convect, area_convect = convect(mld_ws, area_ws, total_mld_ws)

percent_yrs_convect = (ma.count(year_convect))/float(len(year_convect))*100

avg_area_convect = ma.mean(area_convect)/1e5

fig = plt.figure()
plt.plot(total_mld_ws.data, color = 'k', lw = 1.5)
plt.plot( mld_convect, color = '#368BC1', lw = 1.5)
plt.xlabel('Time (years)')
plt.ylabel('meters')
plt.title('(a) Mixed Layer Depth Timeseries', fontsize = 11)
plt.ylim([0,800])
plt.gca().invert_yaxis()

fig = plt.figure()
ax1 = plt.subplot(1,1,1)
t=np.arange(0,500,1)
plt.plot(t, total_mld_ws.data, color = 'k', lw = 1.5)
collection = collections.BrokenBarHCollection.span_where(
        t, ymin=0, ymax=900, where=total_mld_ws.data > 400, facecolor='green', alpha=0.5)
ax1.add_collection(collection)
plt.gca().invert_yaxis()
plt.show()
