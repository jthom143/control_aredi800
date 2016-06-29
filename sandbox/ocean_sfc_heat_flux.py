"""
    Script to calculate surface ocean heat fluxes
    """

import iris
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np
import glob
from netcdftime import utime
from datetime import datetime
import time
import cartopy.crs as ccrs
from netCDF4 import Dataset

plt.ion()

# Today's Date:
timestr = time.strftime('%m%d%Y')

name = '800'

# Import geolat_t
fh = Dataset('/datascope/gnana_esms/jthom143/newCO2_control/ocean_hgrid.nc', mode = 'r')
dy = fh.variables['dy'][:]
dy = dy[:,40]

dy_new = np.zeros(80)

i = 0
for x in range(0,80):
    dy_new[x] = dy[i]+dy[i+1]
    i = i+2


# Load Temp Data
PATHS={'aredi_800':'~/control_aredi800/data/newCO2_control_800/'}



print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    sw_flx = iris.load_cube(PATH+'swflx.nc')
    lw_flx = iris.load_cube(PATH+'lw_heat.nc')
    sens_flx = iris.load_cube(PATH+'sens_heat.nc')
    evap_flx = iris.load_cube(PATH+'evap_heat.nc')
    area = iris.load_cube(PATH+'area_t.nc')


# Zonally Average
sw_flx_zonalavg = sw_flx.collapsed('longitude', iris.analysis.MEAN)
lw_flx_zonalavg = lw_flx.collapsed('longitude', iris.analysis.MEAN)
sens_flx_zonalavg = sens_flx.collapsed('longitude', iris.analysis.MEAN)
evap_flx_zonalavg = evap_flx.collapsed('longitude', iris.analysis.MEAN)

# Time Mean
sw_flx_mean = sw_flx_zonalavg.collapsed('time', iris.analysis.MEAN)
lw_flx_mean = lw_flx_zonalavg.collapsed('time', iris.analysis.MEAN)
sens_flx_mean = sens_flx_zonalavg.collapsed('time', iris.analysis.MEAN)
evap_flx_mean = evap_flx_zonalavg.collapsed('time', iris.analysis.MEAN)
sum = sw_flx_mean + lw_flx_mean + sens_flx_mean + evap_flx_mean

# Zonally Integrate
sw_flx_W = sw_flx*area
lw_flx_W = lw_flx*area
sens_flx_W = sens_flx*area
evap_flx_W = evap_flx*area

sw_flx_integrate = sw_flx_W.collapsed('longitude', iris.analysis.SUM)
lw_flx_integrate = lw_flx_W.collapsed('longitude', iris.analysis.SUM)
sens_flx_integrate = sens_flx_W.collapsed('longitude', iris.analysis.SUM)
evap_flx_integrate = evap_flx_W.collapsed('longitude', iris.analysis.SUM)

sw_flx_integrate = (sw_flx_integrate/dy_new)*111131.7       # W/deg. latitude
lw_flx_integrate = (lw_flx_integrate/dy_new)*111131.7       # W/deg. latitude
sens_flx_integrate = (sens_flx_integrate/dy_new)*111131.7   # W/deg. latitude
evap_flx_integrate = (evap_flx_integrate/dy_new)*111131.7   # W/deg. latitude


# Plot Flux with Latitude for Convective and Non-Convective years
lats = sw_flx_zonalavg.coord('latitude').points

fig1 = plt.figure(figsize = (7.5, 7.5))
ax2 = fig1.add_axes([0.15,0.55,0.7,0.35])
plt.title('(a) Climatology', fontsize = 12)
plt.plot(lats, sw_flx_mean.data, color = '#C11B17', lw = 1.5)
plt.plot(lats, lw_flx_mean.data, color = '#41A317', lw = 1.5)
plt.plot(lats, sens_flx_mean.data, color = '#F87217', lw = 1.5)
plt.plot(lats, evap_flx_mean.data, color = '#306EFF', lw = 1.5)
#plt.plot(lats, sum.data, color = 'k', lw = 1.5)
plt.axhline(0, color = 'k', ls = '--')
plt.xlim([-80, 0])
plt.text(-20, 200, 'Shortwave', color = '#C11B17', fontsize = 9)
plt.text(-20, -90, 'Longwave', color = '#41A317', fontsize = 9)
plt.text(-20, -40, 'Sensible', color = '#F87217', fontsize = 9)
plt.text(-20, -175, 'Latent', color = '#306EFF', fontsize = 9)
plt.text(-75, 200, '* Positive Heats Ocean', color = 'k', fontsize = 10)
plt.ylabel('Ocean Heat Flux [W m$^{-2}$]', fontsize = 12)


ax1 = fig1.add_axes([0.15,0.1,0.7,0.35])
plt.title('(b) Convective - Non-convective', fontsize = 12)
plt.plot(lats, (sw_flx_integrate[237,:].data-sw_flx_integrate[203,:].data)/1e12, color = '#C11B17', lw = 1.5)
plt.plot(lats, (lw_flx_integrate[237,:].data-lw_flx_integrate[203,:].data)/1e12, color = '#41A317', lw = 1.5)
plt.plot(lats, (sens_flx_integrate[237,:].data-sens_flx_integrate[203,:].data)/1e12, color = '#F87217', lw = 1.5)
plt.plot(lats, (evap_flx_integrate[237,:].data-evap_flx_integrate[203,:].data)/1e12, color = '#306EFF', lw = 1.5)
plt.axhline(0, color = 'k', ls = '--')
plt.xlabel('Latitude', fontsize = 12)
plt.ylabel(' Heat Flux [10$^{12}$ W deg$^{-1}$]', fontsize = 12)

sum_int = (sw_flx_integrate[237,:].data-sw_flx_integrate[203,:].data) + (lw_flx_integrate[237,:].data-lw_flx_integrate[203,:].data) + (sens_flx_integrate[237,:].data-sens_flx_integrate[203,:].data) + (evap_flx_integrate[237,:].data-evap_flx_integrate[203,:].data)

fig2 = plt.figure(figsize = (7.5, 7.5))
ax1 = fig2.add_axes([0.15,0.1,0.7,0.35])
plt.title('Sum of Surface Ocean Heat Fluxes', fontsize = 12)
plt.plot(lats, sum_int.data/1e12, color = 'k', lw = 1.5)
plt.xlabel('Latitude', fontsize = 12)
plt.ylabel(' Heat Flux [10$^{12}$ W deg$^{-1}$]', fontsize = 12)
plt.text(-75, 16, '* Positive Heats Ocean', color = 'k', fontsize = 10)
plt.axhline(0, color = 'k', ls = '--')


# Total Heat
total_heat = sw_flx_integrate + lw_flx_integrate + sens_flx_integrate + evap_flx_integrate
total_heat_mean = total_heat.collapsed('time', iris.analysis.MEAN)
total_heat_anomaly = total_heat - total_heat_mean

fig3 = plt.figure(figsize = (7.5, 7.5))
ax2 = fig3.add_axes([0.15,0.55,0.7,0.35])
plt.title('(a) Total Heat Flux', fontsize = 12)
plt.plot(lats, total_heat_anomaly[237].data/1e12, color = 'r', lw = 1.5)
plt.plot(lats, total_heat_anomaly[203].data/1e12, color = 'b', lw = 1.5)
plt.axhline(0, color = 'k', ls = '--')
plt.ylabel(' Heat Flux [10$^{12}$ W deg$^{-1}$]', fontsize = 12)
plt.text(-52, 6, 'Convective', color = 'r', fontsize = 10)
plt.text(-35, -12, 'Non-convective', color = 'b', fontsize = 10)


ax1 = fig3.add_axes([0.15,0.1,0.7,0.35])
plt.title('(b) Convective - Non-convective', fontsize = 12)
plt.plot(lats, sum_int.data/1e12, color = 'k', lw = 1.5)
plt.xlabel('Latitude', fontsize = 12)
plt.ylabel(' Heat Flux [10$^{12}$ W deg$^{-1}$]', fontsize = 12)
plt.axhline(0, color = 'k', ls = '--')

# Individual Fluxes
sw_mean = sw_flx_integrate.collapsed('time', iris.analysis.MEAN)
lw_mean = lw_flx_integrate.collapsed('time', iris.analysis.MEAN)
sens_mean = sens_flx_integrate.collapsed('time', iris.analysis.MEAN)
evap_mean = evap_flx_integrate.collapsed('time', iris.analysis.MEAN)

fig4 = plt.figure(figsize = (7.5, 7.5))
ax2 = fig4.add_axes([0.15,0.79,0.7,0.16])
plt.plot(lats, (sw_flx_integrate[237,:].data - sw_mean.data)/1e12, color = '#C11B17', lw = 1.5)
plt.plot(lats, (sw_flx_integrate[203,:].data - sw_mean.data)/1e12, color = '#C11B17', lw = 1.5, ls = '--')
plt.axhline(0, color = 'k', ls = '--')
plt.title('(a) Shortwave', fontsize = 12)


ax2 = fig4.add_axes([0.15,0.56,0.7,0.16])
plt.plot(lats, (lw_flx_integrate[237,:].data - lw_mean.data)/1e12, color = '#41A317', lw = 1.5)
plt.plot(lats, (lw_flx_integrate[203,:].data - lw_mean.data)/1e12, color = '#41A317', lw = 1.5, ls = '--')
plt.axhline(0, color = 'k', ls = '--')
plt.title('(b) Longwave', fontsize = 12)


ax2 = fig4.add_axes([0.15,0.33,0.7,0.16])
plt.plot(lats, (sens_flx_integrate[237,:].data - sens_mean.data)/1e12, color = '#F87217', lw = 1.5)
plt.plot(lats, (sens_flx_integrate[203,:].data - sens_mean.data)/1e12, color = '#F87217', lw = 1.5, ls = '--')
plt.axhline(0, color = 'k', ls = '--')
plt.title('(c) Sensible', fontsize = 12)
plt.ylabel('                            Heat Flux [10$^{12}$ W deg$^{-1}$]', fontsize = 12)

ax2 = fig4.add_axes([0.15,0.1,0.7,0.16])
plt.plot(lats, (evap_flx_integrate[237,:].data - evap_mean.data)/1e12, color = '#306EFF', lw = 1.5)
plt.plot(lats, (evap_flx_integrate[203,:].data - evap_mean.data)/1e12, color = '#306EFF', lw = 1.5, ls = '--')
plt.axhline(0, color = 'k', ls = '--')
plt.title('(d) Latent', fontsize = 12)
plt.xlabel('Latitude', fontsize = 12)



## Plot convective - non-convective (min-max) for each period
f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
datas = (total_heat[136] - total_heat[99])/1e12
sum = datas.collapsed('latitude', iris.analysis.SUM)
ax1.plot(lats,datas.data , color = 'k')
ax1.set_title('Period 1', fontsize = 10)
ax1.axhline(0, color = 'k', ls = '--')
ax1.set_xlim([-80,80])
ax1.text(-75, 17, 'Sum = %d W' % sum.data, color = 'k', fontsize = 10)

datas = (total_heat[237] - total_heat[203])/1e12
sum2 = datas.collapsed('latitude', iris.analysis.SUM)
ax2.plot(lats, datas.data, color = 'k')
ax2.set_title('Period 2', fontsize = 10)
ax2.axhline(0, color = 'k', ls = '--')
ax2.set_xlim([-80,80])
ax2.text(-75, 17, 'Sum = %d W' % sum2.data, color = 'k', fontsize = 10)

datas = (total_heat[297] - total_heat[277])/1e12
sum3 = datas.collapsed('latitude', iris.analysis.SUM)
ax3.plot(lats,datas.data, color = 'k')
ax3.set_title('Period 3', fontsize = 10)
ax3.axhline(0, color = 'k', ls = '--')
ax3.set_xlim([-80,80])
ax3.set_ylabel('Heat Flux [10$^{12}$ W deg$^{-1}$]', fontsize = 10)
ax3.text(-75, 28, 'Sum = %d W' % sum3.data, color = 'k', fontsize = 10)

datas = (total_heat[352] - total_heat[330])/1e12
sum4 = datas.collapsed('latitude', iris.analysis.SUM)
ax4.plot(lats,datas.data, color = 'k')
ax4.set_title('Period 4', fontsize = 10)
ax4.axhline(0, color = 'k', ls = '--')
ax4.set_xlim([-80,80])
ax4.text(-75, 28, 'Sum = %d W' % sum4.data, color = 'k', fontsize = 10)

datas = (total_heat[405] - total_heat[380])/1e12
sum5 = datas.collapsed('latitude', iris.analysis.SUM)
ax5.plot(lats, datas.data, color = 'k')
ax5.set_title('Period 5', fontsize = 10)
ax5.axhline(0, color = 'k', ls = '--')
ax5.set_xlim([-80,80])
ax5.text(-75, 28, 'Sum = %d W' % sum5.data, color = 'k', fontsize = 10)

datas = (total_heat[499] - total_heat[462])/1e12
sum6 = datas.collapsed('latitude', iris.analysis.SUM)
ax6.plot(lats, datas.data, color = 'k')
ax6.set_title('Period 6', fontsize = 10)
ax6.axhline(0, color = 'k', ls = '--')
ax6.set_xlim([-80,80])
ax6.text(-75, 28, 'Sum = %d W' % sum6.data, color = 'k', fontsize = 10)
plt.savefig('notes/figures/total_heat_flux_anomaly.png')


## Plot convective - non-convective (min-max) for each period SOUTHERN HEMISPHERE
f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')

ax1.plot(lats,(total_heat[136].data - total_heat[99].data)/1e12, color = 'k')
ax1.set_title('Period 1', fontsize = 10)
ax1.axhline(0, color = 'k', ls = '--')
ax1.set_xlim([-80,0])

ax2.plot(lats,(total_heat[237].data - total_heat[203].data)/1e12, color = 'k')
ax2.set_title('Period 2', fontsize = 10)
ax2.axhline(0, color = 'k', ls = '--')
ax2.set_xlim([-80,0])


ax3.plot(lats,(total_heat[297].data - total_heat[277].data)/1e12, color = 'k')
ax3.set_title('Period 3', fontsize = 10)
ax3.axhline(0, color = 'k', ls = '--')
ax3.set_xlim([-80,0])
ax3.set_ylabel('Heat Flux [10$^{12}$ W deg$^{-1}$]', fontsize = 10)


ax4.plot(lats,(total_heat[352].data - total_heat[330].data)/1e12, color = 'k')
ax4.set_title('Period 4', fontsize = 10)
ax4.axhline(0, color = 'k', ls = '--')
ax4.set_xlim([-80,0])


ax5.plot(lats,(total_heat[405].data - total_heat[380].data)/1e12, color = 'k')
ax5.set_title('Period 5', fontsize = 10)
ax5.axhline(0, color = 'k', ls = '--')
ax5.set_xlim([-80,0])


ax6.plot(lats,(total_heat[499].data - total_heat[462].data)/1e12, color = 'k')
ax6.set_title('Period 6', fontsize = 10)
ax6.axhline(0, color = 'k', ls = '--')
ax6.set_xlim([-80,0])


## Plot convective - non-convective (min-max) for each period SHORTWAVE
f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
ax1.plot(lats,(total_heat[136].data - total_heat[99].data)/1e12, color = 'k')
ax1.plot(lats,(sw_flx_integrate[136].data - sw_flx_integrate[99].data)/1e12 , color = '#C11B17')
ax1.plot(lats,(lw_flx_integrate[136].data - lw_flx_integrate[99].data)/1e12 , color = '#41A317')
ax1.plot(lats,(sens_flx_integrate[136].data - sens_flx_integrate[99].data)/1e12 , color = '#F87217')
ax1.plot(lats,(evap_flx_integrate[136].data - evap_flx_integrate[99].data)/1e12 , color = '#306EFF')
ax1.set_title('Period 1', fontsize = 10)
ax1.axhline(0, color = 'k', ls = '--')
ax1.set_xlim([-80,80])

ax2.plot(lats,(total_heat[237].data - total_heat[203].data)/1e12, color = 'k')
ax2.plot(lats, (sw_flx_integrate[237].data - sw_flx_integrate[203].data)/1e12, color = '#C11B17')
ax2.plot(lats, (lw_flx_integrate[237].data - lw_flx_integrate[203].data)/1e12, color = '#41A317')
ax2.plot(lats, (sens_flx_integrate[237].data - sens_flx_integrate[203].data)/1e12, color = '#F87217')
ax2.plot(lats, (evap_flx_integrate[237].data - evap_flx_integrate[203].data)/1e12, color = '#306EFF')
ax2.set_title('Period 2', fontsize = 10)
ax2.axhline(0, color = 'k', ls = '--')
ax2.set_xlim([-80,80])

ax3.plot(lats,(total_heat[297].data - total_heat[277].data)/1e12, color = 'k')
ax3.plot(lats,(sw_flx_integrate[297].data - sw_flx_integrate[277].data)/1e12, color = '#C11B17')
ax3.plot(lats,(lw_flx_integrate[297].data - lw_flx_integrate[277].data)/1e12, color = '#41A317')
ax3.plot(lats,(sens_flx_integrate[297].data - sens_flx_integrate[277].data)/1e12, color = '#F87217')
ax3.plot(lats,(evap_flx_integrate[297].data - evap_flx_integrate[277].data)/1e12, color = '#306EFF')
ax3.set_title('Period 3', fontsize = 10)
ax3.axhline(0, color = 'k', ls = '--')
ax3.set_xlim([-80,80])
ax3.set_ylabel('Heat Flux [10$^{12}$ W deg$^{-1}$]', fontsize = 10)

ax4.plot(lats,(total_heat[352].data - total_heat[330].data)/1e12, color = 'k')
ax4.plot(lats,(sw_flx_integrate[352].data - sw_flx_integrate[330].data)/1e12, color = '#C11B17')
ax4.plot(lats,(lw_flx_integrate[352].data - lw_flx_integrate[330].data)/1e12, color = '#41A317')
ax4.plot(lats,(sens_flx_integrate[352].data - sens_flx_integrate[330].data)/1e12, color = '#F87217')
ax4.plot(lats,(evap_flx_integrate[352].data - evap_flx_integrate[330].data)/1e12, color = '#306EFF')
ax4.set_title('Period 4', fontsize = 10)
ax4.axhline(0, color = 'k', ls = '--')
ax4.set_xlim([-80,80])

ax5.plot(lats,(total_heat[405].data - total_heat[380].data)/1e12, color = 'k')
ax5.plot(lats, (sw_flx_integrate[405].data - sw_flx_integrate[380].data)/1e12, color = '#C11B17')
ax5.plot(lats, (lw_flx_integrate[405].data - lw_flx_integrate[380].data)/1e12, color = '#41A317')
ax5.plot(lats, (sens_flx_integrate[405].data - sens_flx_integrate[380].data)/1e12, color = '#F87217')
ax5.plot(lats, (evap_flx_integrate[405].data - evap_flx_integrate[380].data)/1e12, color = '#306EFF')
ax5.set_title('Period 5', fontsize = 10)
ax5.axhline(0, color = 'k', ls = '--')
ax5.set_xlim([-80,80])

ax6.plot(lats,(total_heat[499].data - total_heat[462].data)/1e12, color = 'k')
ax6.plot(lats, (sw_flx_integrate[499].data - sw_flx_integrate[462].data)/1e12, color = '#C11B17')
ax6.plot(lats, (lw_flx_integrate[499].data - lw_flx_integrate[462].data)/1e12, color = '#41A317')
ax6.plot(lats, (sens_flx_integrate[499].data - sens_flx_integrate[462].data)/1e12, color = '#F87217')
ax6.plot(lats, (evap_flx_integrate[499].data - evap_flx_integrate[462].data)/1e12, color = '#306EFF')
ax6.set_title('Period 6', fontsize = 10)
ax6.axhline(0, color = 'k', ls = '--')
ax6.set_xlim([-80,80])
plt.savefig('notes/figures/heat_flux_components_anomaly.png')



###
f, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex='col', sharey='row')
ax1.plot(lats,(total_heat[136].data - total_heat[99].data)/1e12, color = 'k')
ax1.plot(lats,((sw_flx_integrate[136].data - sw_flx_integrate[99].data)/1e12)+((evap_flx_integrate[136].data - evap_flx_integrate[99].data)/1e12) , color = '#C11B17')
ax1.set_title('Period 1', fontsize = 10)
ax1.axhline(0, color = 'k', ls = '--')
ax1.set_xlim([-80,80])

ax2.plot(lats,(total_heat[237].data - total_heat[203].data)/1e12, color = 'k')
ax2.plot(lats, ((sw_flx_integrate[237].data - sw_flx_integrate[203].data)/1e12)+((evap_flx_integrate[237].data - evap_flx_integrate[203].data)/1e12), color = '#C11B17')
ax2.set_title('Period 2', fontsize = 10)
ax2.axhline(0, color = 'k', ls = '--')
ax2.set_xlim([-80,80])

ax3.plot(lats,(total_heat[297].data - total_heat[277].data)/1e12, color = 'k')
ax3.plot(lats,((sw_flx_integrate[297].data - sw_flx_integrate[277].data)/1e12)+((evap_flx_integrate[297].data - evap_flx_integrate[277].data)/1e12), color = '#C11B17')
ax3.set_title('Period 3', fontsize = 10)
ax3.axhline(0, color = 'k', ls = '--')
ax3.set_xlim([-80,80])
ax3.set_ylabel('Heat Flux [10$^{12}$ W deg$^{-1}$]', fontsize = 10)

ax4.plot(lats,(total_heat[352].data - total_heat[330].data)/1e12, color = 'k')
ax4.plot(lats,((sw_flx_integrate[352].data - sw_flx_integrate[330].data)/1e12)+((evap_flx_integrate[352].data - evap_flx_integrate[330].data)/1e12), color = '#C11B17')
ax4.set_title('Period 4', fontsize = 10)
ax4.axhline(0, color = 'k', ls = '--')
ax4.set_xlim([-80,80])

ax5.plot(lats,(total_heat[405].data - total_heat[380].data)/1e12, color = 'k')
ax5.plot(lats, ((sw_flx_integrate[405].data - sw_flx_integrate[380].data)/1e12)+((evap_flx_integrate[405].data - evap_flx_integrate[380].data)/1e12), color = '#C11B17')
ax5.set_title('Period 5', fontsize = 10)
ax5.axhline(0, color = 'k', ls = '--')
ax5.set_xlim([-80,80])

ax6.plot(lats,(total_heat[499].data - total_heat[462].data)/1e12, color = 'k')
ax6.plot(lats, ((sw_flx_integrate[499].data - sw_flx_integrate[462].data)/1e12)+((evap_flx_integrate[499].data - evap_flx_integrate[462].data)/1e12), color = '#C11B17')
ax6.set_title('Period 6', fontsize = 10)
ax6.axhline(0, color = 'k', ls = '--')
ax6.set_xlim([-80,80])
plt.savefig('notes/figures/heat_flux_sw+evap_anomaly.png')






