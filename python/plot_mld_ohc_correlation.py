"""
    Script to plot the correlation between convection and heat content using seaborn
"""
import iris
import numpy as np
import seaborn as sns
sns.set(style="darkgrid", color_codes=True)

# Load Data
names = {'aredi_400', 'aredi_800', 'aredi_2400'}

mld = {}
global_heat = {}

for name in names:
    mld[name] = np.load('/RESEARCH/control_aredi800/data/tmp_npy/'+name+'mld.npy')
    global_heat[name] = np.load('/RESEARCH/control_aredi800/data/tmp_npy/'+name+'global_heat.npy')

x = mld['aredi_800']
y = global_heat['aredi_800']
g = sns.jointplot(x, y, kind="reg")
