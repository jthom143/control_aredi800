import iris
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
import glob
from netcdftime import utime
from datetime import datetime
import time
import cartopy.crs as ccrs
import scipy.stats as stats

start_time = time.time()

def calc_fluxes(fluxes_ann, fluxes, area, domain):

    
    # Area of Interest
    constraint = iris.Constraint(longitude=lambda x: domain[0] < x < domain[1], latitude=lambda y: domain[2] < y < domain[3], depth=lambda z: domain[4] < z < domain[5])
    constraint_area = iris.Constraint(longitude=lambda x: domain[0] < x < domain[1], latitude=lambda y: domain[2] < y < domain[3])
    t = 0
    
    # Calculate Average Fluxes in area of interest
        
    for name, var in fluxes.iteritems():
        if name == 'temp_xflux_adv':
            fluxes_ann[name] = fluxes_ann[name].extract(constraint)
            x_flux = fluxes_ann[name][t,:,:,:].collapsed(['depth', 'latitude'], iris.analysis.MEAN)
            fluxes_ann[name] = x_flux[0] - x_flux[-1]
        elif name == 'temp_yflux_adv':
            fluxes_ann[name] = fluxes_ann[name].extract(constraint)
            y_flux = fluxes_ann[name][t,:,:,:].collapsed(['depth', 'longitude'], iris.analysis.MEAN)
            fluxes_ann[name] = y_flux[0] - y_flux[-1]             
        elif name =='temp_zflux_adv':
            fluxes_ann[name] = fluxes_ann[name].extract(constraint)
            z_flux = fluxes_ann[name][t,:,:,:].collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
            fluxes_ann[name] = z_flux[-1] - z_flux[-0]             
        else:
            fluxes_ann[name] = fluxes_ann[name].extract(constraint)
            fluxes_ann[name] = fluxes_ann[name][t,:,:,:].collapsed(['depth', 'latitude', 'longitude'], iris.analysis.SUM)                                         
    area_domain = area.extract(constraint_area)
    area_domain = area_domain.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    return fluxes_ann, area_domain



CLIM_PATH = '/datascope/pradal-esm/newCO2_control_800/history/'
clim_files = sorted(glob.glob(CLIM_PATH+'*ocean_month.nc'))

fluxes= {'temp_xflux_adv':'cp*rho_dxt*dyt_u_temp',
         'temp_yflux_adv':'cp*rho*dzt*dxt*v*temp',
         'temp_zflux_adv':'cp*rho*dxt*dyt*wt*temp',
         'temp_tendency':'time tendency for tracer Potential temperature',
         'neutral_temp':'rho*dzt*cp*explicit neutral tendency (heating)',
         'temp_submeso':'rho*dzt*cp*submesoscale tendency (heating)',
         'temp_nonlocal_KPP': 'cp*rho*dzt*nonlocal tendency from KPP',
         'temp_vdiffuse_impl':'implicit vert diffusion of heat',
         'sw_heat':'downwelling_shortwave_flux_in_sea_water'
        }


raw_flux = {}

for name, var in fluxes.iteritems():
    raw_flux[name] = iris.load_cube('~/control_aredi800/data/'+name+'.nc')
    if name =='temp_zflux_adv':
        raw_flux[name].coord('ucell pstar').standard_name = 'depth'
    else:
        raw_flux[name].coord('tcell pstar').standard_name = 'depth'
    

area = iris.load_cube(clim_files[0], 'tracer cell area')
area.remove_coord(area.aux_coords[0])
area.remove_coord(area.aux_coords[0])
hu = iris.load_cube(clim_files[0], 'ocean u-cell thickness')[0,:]
hu.remove_coord(hu.aux_coords[0])
hu.remove_coord(hu.aux_coords[0])
hu.remove_coord(hu.aux_coords[0])
ht = iris.load_cube(clim_files[0], 'cell_thickness')[0,:]
ht.remove_coord(ht.aux_coords[0])
ht.remove_coord(ht.aux_coords[0])
ht.remove_coord(ht.aux_coords[0])

domain =  [-60, -40,-70, -50, 100, 500 ]
fluxes_ann, area_domain = calc_fluxes(raw_flux, fluxes, area, domain)

temp_tendency = fluxes_ann['temp_tendency'].data
advective_flux = (fluxes_ann['temp_xflux_adv'].data + fluxes_ann['temp_yflux_adv'].data + fluxes_ann['temp_zflux_adv'].data)/area_domain.data
eddy_flux = fluxes_ann['temp_submeso'].data + fluxes_ann['neutral_temp'].data
turb_flux = fluxes_ann['temp_nonlocal_KPP'].data + fluxes_ann['temp_vdiffuse_impl'].data + fluxes_ann['sw_heat'].data
SUM = advective_flux+eddy_flux+turb_flux

print 'Temperature Tendency: %f' % temp_tendency
print 'Advective Flux: %f' % advective_flux
print 'Eddy Flux: %f ' % eddy_flux
print 'Turbulent Flux: %f' % turb_flux
print 'SUM = %f' %SUM


print '   '
print '   '
print("--- %s seconds ---" % (time.time() - start_time))
