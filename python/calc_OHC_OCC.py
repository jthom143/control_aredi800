"""
Script to calculate ocean heat and carbon content
"""

import iris
import matplotlib.pyplot as plt
import numpy as np

def calc_ohc(temp, rhodzt, area):

    # Convert Units and multiply by rho cp to change units from C to J/m^2
    temp = temp + 274.15
    temp = temp * rhodzt
    heat = temp * 4000
    heat.units = 'J m^-2'

    # Vertically sum
    heat_sumz = heat.collapsed('tcell pstar', iris.analysis.SUM)

    # Sum over globe
    heat_sumz_weighted = heat_sumz*area
    heat_global = heat_sumz_weighted.collapsed(['longitude', 'latitude'], iris.analysis.SUM)

    return heat, heat_sumz, heat_global

def calc_occ(dic, rhodzt, area):

    # Convert Units and multiply by rho cp to change units from g to J/m^2
    dic = dic * rhodzt
    dic = dic * 12 
    dic.units = 'g m^-2'

    # Vertically Sum
    dic_sumz = dic.collapsed('tcell pstar', iris.analysis.SUM)

    # Sum over globe
    dic_sumz_weighted = dic_sumz*area
    dic_global = dic_sumz_weighted.collapsed(['longitude', 'latitude'], iris.analysis.SUM)

    return dic, dic_sumz, dic_global



