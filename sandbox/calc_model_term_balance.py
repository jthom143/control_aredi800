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


def calc_fluxes(fluxes_ann, fluxes, area, ht):
    
    # Grid Box of Interest
    i = 80 # Longitude
    j = 11 # Latitude
    k = 7  # Depth
    
    # Print Sample Box Coordinates:
    print fluxes_ann['temp_xflux_adv'].coord('longitude').points[i]
    print fluxes_ann['temp_yflux_adv'].coord('latitude').points[j]
    print fluxes_ann['temp_zflux_adv'].coord('ucell pstar').points[k]

    for name, var in fluxes.iteritems():
        if name == 'temp_xflux_adv':
            fluxes_ann[name] = fluxes_ann[name][:,k,j,i-1] - fluxes_ann[name][:,k,j,i]
        elif name == 'temp_yflux_adv':
            fluxes_ann[name] = fluxes_ann[name][:,k,j-1,i] - fluxes_ann[name][:,k,j,i]             
        elif name =='temp_zflux_adv':
            fluxes_ann[name] = fluxes_ann[name][:,k,j,i] - fluxes_ann[name][:,k-1,j,i]             
        else:
            fluxes_ann[name] = fluxes_ann[name][:,k,j,i]                                         

    return fluxes_ann



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
    print name



area = iris.load_cube(clim_files[0], 'tracer cell area')
hu = iris.load_cube(clim_files[0], 'ocean u-cell thickness')[0,:]
hu.remove_coord(hu.aux_coords[0])
hu.remove_coord(hu.aux_coords[0])
hu.remove_coord(hu.aux_coords[0])
ht = iris.load_cube(clim_files[0], 'cell_thickness')[0,:]
ht.remove_coord(ht.aux_coords[0])
ht.remove_coord(ht.aux_coords[0])
ht.remove_coord(ht.aux_coords[0])


fluxes_ann = calc_fluxes(raw_flux, fluxes, area, ht)

time = fluxes_ann['temp_tendency'].coord('time').points

fig1 = plt.figure()
plt.plot(fluxes_ann['temp_tendency'].data, color = 'k', lw = 1.5 ,label = 'Temperature Tendency')
plt.plot((fluxes_ann['temp_submeso'].data+fluxes_ann['neutral_temp'].data), color = 'g', label = 'eddy flux')
plt.plot((fluxes_ann['temp_vdiffuse_impl'].data + fluxes_ann['temp_nonlocal_KPP'].data + fluxes_ann['sw_heat'].data), color = 'r', label = 'turbulence flux')
plt.plot((fluxes_ann['temp_xflux_adv'].data+fluxes_ann['temp_yflux_adv'].data + fluxes_ann['temp_zflux_adv'].data)/area[11, 80].data, color = 'b', label = 'advective flux')
plt.plot((fluxes_ann['temp_xflux_adv'].data+fluxes_ann['temp_yflux_adv'].data + fluxes_ann['temp_zflux_adv'].data)/area[11,80].data +fluxes_ann['temp_submeso'].data+fluxes_ann['neutral_temp'].data+fluxes_ann['temp_vdiffuse_impl'].data + fluxes_ann['temp_nonlocal_KPP'].data + fluxes_ann['sw_heat'].data, ls = '--', color = 'k', label = 'Sum')
plt.xlim([0,100])
plt.xlabel('Time (years)')
plt.ylabel("W m$\mathrm{^{-2}}$")
plt.legend()
plt.show()
