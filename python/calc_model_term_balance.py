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

#data_PATH = '~/control_aredi800/data/'
data_PATH = '~/Desktop/data/'

def calc_fluxes(fluxes_ann, fluxes, area, domain):

    
    # Area of Interest:
    constraint = iris.Constraint(longitude=lambda x: domain[0] < x < domain[1], latitude=lambda y: domain[2] < y < domain[3], depth=lambda z: domain[4] < z < domain[5])
    constraint_area = iris.Constraint(longitude=lambda x: domain[0] < x < domain[1], latitude=lambda y: domain[2] < y < domain[3])

    # Time of Interest:
    t = 0    

    # Calculate Average Fluxes in area of interest
    '''    
    for name, var in fluxes.iteritems():
        if name == 'temp_xflux_adv':
            fluxes_ann[name] = fluxes_ann[name].extract(constraint)
            x_flux = fluxes_ann[name][t,:,:,:].collapsed(['depth', 'latitude'], iris.analysis.SUM)
            fluxes_ann[name] = x_flux[0] - x_flux[-1]
        elif name == 'temp_yflux_adv':
            fluxes_ann[name] = fluxes_ann[name].extract(constraint)
            y_flux = fluxes_ann[name][t,:,:,:].collapsed(['depth', 'longitude'], iris.analysis.SUM)
            fluxes_ann[name] = y_flux[0] - y_flux[-1]             
        elif name =='temp_zflux_adv':
            fluxes_ann[name] = fluxes_ann[name].extract(constraint)
            z_flux = fluxes_ann[name][t,:,:,:].collapsed(['latitude', 'longitude'], iris.analysis.SUM)
            fluxes_ann[name] = z_flux[-1] - z_flux[-0]             
        else:
            fluxes_ann[name] = fluxes_ann[name].extract(constraint)
            fluxes_ann[name] = fluxes_ann[name][t,1:,1:,1:].collapsed(['depth', 'latitude', 'longitude'], iris.analysis.SUM)
            
    area_domain = area.extract(constraint_area)
    area_domain = area_domain[:,:].collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    return fluxes_ann, area_domain
    '''
    for name, var in fluxes.iteritems():
        fluxes_ann[name] = fluxes_ann[name][t,:,:,:].extract(constraint)

        fluxes_domain= {'temp_xflux_adv':np.full_like(fluxes_ann['temp_tendency'].data, 1),
                        'temp_yflux_adv':np.full_like(fluxes_ann['temp_tendency'].data, 1),
                        'temp_zflux_adv':np.full_like(fluxes_ann['temp_tendency'].data, 1),
                        'temp_tendency':np.full_like(fluxes_ann['temp_tendency'].data, 1),
                        'neutral_temp':np.full_like(fluxes_ann['temp_tendency'].data, 1),
                        'temp_submeso':np.full_like(fluxes_ann['temp_tendency'].data, 1),
                        'temp_nonlocal_KPP':np.full_like(fluxes_ann['temp_tendency'].data, 1),
                        'temp_vdiffuse_impl':np.full_like(fluxes_ann['temp_tendency'].data, 1),
                        'sw_heat':np.full_like(fluxes_ann['temp_tendency'].data, 1)
                        }
    
    for name, var in fluxes.iteritems():
        print name

        for i in range(0,len(fluxes_ann['temp_vdiffuse_impl'].coord('longitude').points)-1):
            for j in range(0,len(fluxes_ann['temp_vdiffuse_impl'].coord('latitude').points)-1):
                for k in range(0,len(fluxes_ann['temp_vdiffuse_impl'].coord('depth').points)-1):
                    if name =='temp_xflux_adv':
                        fluxes_domain[name][k,j,i] = fluxes_ann[name][k,j,i].data - fluxes_ann[name][k,j,i+1].data
                    elif name == 'temp_yflux_adv':
                        fluxes_domain[name][k,j,i] = fluxes_ann[name][k,j,i].data - fluxes_ann[name][k,j+1,i].data
                    elif name =='temp_zflux_adv':
                        fluxes_domain[name][k,j,i] = fluxes_ann[name][k+1,j,i].data - fluxes_ann[name][k,j,i].data
                    else:
                        fluxes_domain[name][k,j,i] = fluxes_ann[name][k+1, j+1, i+1].data

    

    # Calculate Area of Domain:
    area_domain = area.extract(constraint_area)
    area_domain = area_domain[:,:].collapsed(['longitude', 'latitude'], iris.analysis.SUM)

    return fluxes_domain, area_domain


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
    raw_flux[name] = iris.load_cube(data_PATH+name+'.nc')
    if name =='temp_zflux_adv':
        raw_flux[name].coord('ucell pstar').standard_name = 'depth'
    else:
        raw_flux[name].coord('tcell pstar').standard_name = 'depth'
    

area = iris.load_cube(data_PATH+'area_t.nc')

domain =  [-70, -50,-70, -50, 100, 500 ]
fluxes_ann, area_domain = calc_fluxes(raw_flux, fluxes, area, domain)

'''
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
'''

print '   '
print '   '
print("--- %s seconds ---" % (time.time() - start_time))
