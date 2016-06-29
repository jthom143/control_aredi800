
"""
    Script to compare the 1500-2000m Weddell Sea Temperature for Aredi=800 and GM_min=600 simulations
    """

import iris
import iris.analysis.cartography
import matplotlib.pyplot as plt
import numpy as np
import numpy as np
from scipy.stats import genextreme
import numpy.ma as ma

import sys # access system routines
sys.path.append('/home/jthom143/python_functions')
import colormaps as cmaps


import sys # access system routines
sys.path.append('~/control_aredi800/python')
from calc_region import calc_weddell_sea


## Load Data ###
'''
    Files loaded in this script were created in 'plot_aredi_climatologies.py'
    '''

# Load Temp Data
PATHS={'aredi_800':'~/control_aredi800/data/newCO2_control_800/'}

temp = {}
mld = {}
mld_ws = {}
temp_ws = {}

print 'loading data for:'
for name, PATH in PATHS.iteritems():
    print name
    temp[name] = iris.load_cube(PATH+'temp.nc')
    mld[name] = iris.load_cube(PATH+'mld.nc')
    
    temp_ws[name]= calc_weddell_sea(temp[name])
    mld_ws[name] = calc_weddell_sea(mld[name])
    
    if name == 'aredi_800':
        temp[name] = temp[name][-500:]

    
    mld_ws[name].coord('latitude').guess_bounds()
    mld_ws[name].coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mld_ws[name])
    mld_ws[name] = mld_ws[name].collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    
    temp_ws[name].coord('latitude').guess_bounds()
    temp_ws[name].coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(temp_ws[name])
    temp_ws[name] = temp_ws[name].collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)

mean = np.mean(temp_ws['aredi_800'][:,18].data)
std = np.std(temp_ws['aredi_800'][:,18].data)
convect = mean-std


# Create Figure
plt.figure()
time = np.arange(0, 500, 1)
plt.plot(time, temp_ws['aredi_800'][:,18].data, color = 'k', lw = 1.5, label = 'Aredi = 800')
plt.xlim([0,500])
plt.xlabel('Time [years]')
plt.ylabel('Temperature [$^o$C]')
plt.title('1400m Weddell Sea Temperature')
plt.legend(fontsize = 10)

# PDFS

plt.figure()

plt.hist(temp_ws['aredi_800'][:,18].data, 50, color = 'g', normed=True)
plt.title('Aredi = 800 m$^2$s$^{-1}$', fontsize = 11)

x = np.linspace(0,2.5,100)
param = genextreme.fit(temp_ws['aredi_800'][:,18].data) # distribution fitting
pdf_fitted_800 = genextreme.pdf(x,c = param[0], loc=param[1],scale=param[2])
plt.plot(x,pdf_fitted_800,'k-', lw = 1.5)

plt.axvline(convect, color = 'k', ls = '--', lw = 1.5)


# Timeseries highlighting convective years
plt.figure(figsize=(10,4))
time = np.arange(0, 500, 1)
plt.plot(time, temp_ws['aredi_800'][:,18].data, color = 'k', lw = 1.5, label = 'Aredi = 800')
plt.xlim([0,500])
plt.xlabel('Time [years]')
plt.ylabel('Temperature [$^o$C]')
plt.title('1400m Weddell Sea Temperature')

convective_years = ma.masked_greater(temp_ws['aredi_800'][:,18].data, convect)
plt.plot(time, convective_years, color = 'b', lw = 1.5, label = 'Aredi = 800')

no = ma.count(convective_years)
no2 = ma.count_masked(convective_years)
plt.text(320, 1.7, '# of convecting years = %d' % no, fontsize = 10, color = 'b')
plt.text(320, 1.6, '# of non-convecting Years = %d' % no2, fontsize = 10)



# Histogram on convection
plt.figure()

plt.hist(convective_years.compressed(), 20, color = 'b', normed=True)
plt.hist(temp_ws['aredi_800'][:,18].data, 50, color = 'g')
plt.hist(convective_years.compressed(), 20, color = 'b')
plt.title('convective years', fontsize = 11)
plt.show()

