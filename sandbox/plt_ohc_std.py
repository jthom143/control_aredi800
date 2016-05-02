import iris
import iris.analysis.cartography
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np


from calc_OHC_OCC import calc_ohc
from calc_region import calc_six_regions


# Load Temp Data
PATHS={'aredi_800':'~/control_aredi800/data/newCO2_control_800/'}
       
temp = {}
rhodzt = {}
area = {}
heat = {}
heat_sumz = {}
heat_global = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
        print name
        temp[name] = iris.load_cube(PATH+'temp.nc')
        rhodzt[name] = iris.load_cube(PATH+'rho_dzt.nc')
        area[name] = iris.load_cube(PATH+'area_t.nc')
	heat[name], heat_sumz[name], heat_global[name] = calc_ohc(temp[name], rhodzt[name], area[name])

heat_variance = heat_sumz['aredi_800'].collapsed('time', iris.analysis.VARIANCE)
