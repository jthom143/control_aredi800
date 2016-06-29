## Script to create data used in generating figures for FESD 2016 presentation

import iris
import numpy as np
import sys

# Import my functions
sys.path.append('~/control_aredi800/python')
from calc_region import calc_weddell_sea
from calc_OHC_OCC import calc_occ, calc_ohc

# Load in annually averaged netcdf files from '/home/jthom143/control_aredi800/data/newCO2_control_800/'
PATH = '~/control_aredi800/data/newCO2_control_800/'

temp = iris.load_cube(PATH+'temp.nc')     # Ocean Temperature
rhodzt = iris.load_cube(PATH+'rho_dzt.nc')  # rho*dz/dt
area = iris.load_cube(PATH+'area_t.nc')     # Tracer Cell Area
mld = iris.load_cube(PATH+'mld.nc')         # Mixed Layer Depth
dic = iris.load_cube(PATH+'dic.nc')         # Dissolved Inorganic Carbon

# Trim Ocean Bling Data to match other data lengths
dic = dic[-500:]

# Isolate and average over the the Weddell Sea
temp_ws= calc_weddell_sea(temp)
mld_ws = calc_weddell_sea(mld)

mld_ws.coord('latitude').guess_bounds()
mld_ws.coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(mld_ws)
mld_ws = mld_ws.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    
temp_ws.coord('latitude').guess_bounds()
temp_ws.coord('longitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(temp_ws)
temp_ws = temp_ws.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)

#### Calculate Carbon Content #######################################################################
carbon, carbon_sumz, carbon_global = calc_occ(dic, rhodzt, area)

# Southern Hemisphere Carbon Content:
constraint = iris.Constraint(latitude=lambda y: -90 < y < -50)
carbon_SH_tmp = carbon_sumz.extract(constraint)
area_SH = area.extract(constraint)
b = carbon_SH_tmp*area_SH
carbon_SH = b.collapsed(['longitude', 'latitude'], iris.analysis.SUM)

#########################################################################################################

#### Calculate Heat Content #############################################################################
heat, heat_sumz, heat_global = calc_ohc(temp, rhodzt, area)

# Southern Hemisphere:
constraint = iris.Constraint(latitude=lambda y: -90 < y < -50)
heat_SH_tmp = heat_sumz.extract(constraint)
area_SH = area.extract(constraint)
b = heat_SH_tmp*area_SH
heat_SH = b.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
#########################################################################################################

#### Calculate Greenhouse ###############################################################################
olr = iris.load_cube(PATH+'olr.nc')
lwup_sfc = iris.load_cube(PATH+'lwup_sfc.nc')

G = lwup_sfc - olr
#########################################################################################################


# Save Variables
iris.save(temp_ws, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/temp_ws.nc')
iris.save(mld_ws, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/mld_ws.nc')
iris.save(carbon, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/carbon.nc')
iris.save(heat, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/heat.nc')
iris.save(carbon_sumz, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/carbon_sumz.nc')
iris.save(heat_sumz, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/heat_sumz.nc')
iris.save(G, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/G.nc')
iris.save(carbon_global, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/carbon_global.nc')
iris.save(carbon_SH, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/carbon_SH.nc')
iris.save(heat_global, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/heat_global.nc')
iris.save(heat_SH, '/home/jthom143/control_aredi800/FESD_meeting_2016/data/heat_SH.nc')

