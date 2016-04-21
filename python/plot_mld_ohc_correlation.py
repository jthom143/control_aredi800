"""
    Script to plot the correlation between convection and heat content using seaborn
"""
import iris
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
sns.set(style="white", color_codes=True)


def xcorr(x, y=None, maxlags=None, scale_type='normalize'):
    """Compute the cross correlation of `x` and `y`.

    This function computes the cross correlation of `x` and `y`. It uses the
    equivalence of the cross correlation with the negative convolution computed
    using a FFT to achieve much faster cross correlation than is possible with
    the traditional signal processing definition.

    By default it computes the cross correlation at each of 1 - maxlags to maxlags,
    scaled by the lag 0 cross correlation after mean centering the data.

    Note that it is not necessary for `x` and `y` to be the same size."""






# Load Data
names = {'aredi_400', 'aredi_800', 'aredi_2400'}

mld = {}
global_heat = {}
sh_heat = {}

for name in names:
    mld[name] = np.load('/RESEARCH/control_aredi800/data/tmp_npy/'+name+'_mld.npy')
    global_heat[name] = np.load('/RESEARCH/control_aredi800/data/tmp_npy/'+name+'_global_heat.npy')
    global_heat[name] = global_heat[name]/1e22
    sh_heat[name] = np.load('/RESEARCH/control_aredi800/data/tmp_npy/'+name+'_sh_heat.npy')
    sh_heat[name] = sh_heat[name]/1e22


x = mld['aredi_400']
y = global_heat['aredi_400']
sns.jointplot(x, y, color = 'g', space=0)
#plt.savefig('/RESEARCH/control_aredi800/notes/figures/aredi400_jointplot.png', bbox_inches='tight')

fig = plt.figure()
plt.xcorr(x, y, maxlags = 50, normed = True, color = 'g')
plt.xlabel('Time Lag [years]')
plt.ylabel('Correlation Coefficient')
#plt.savefig('/RESEARCH/control_aredi800/notes/figures/aredi400_xcorr.png', bbox_inches='tight')

plt.show()
