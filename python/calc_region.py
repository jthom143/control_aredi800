"""
Script to calculate variable for specific regions in ocean
"""

def calc_weddell_sea(var):
    import iris
    
    constraint = iris.Constraint(latitude=lambda y: -85 < y < -60, longitude=lambda x: -60 < x < 0)
    var_ws = var.extract(constraint)

    return var_ws


def calc_six_regions(var):
    # Function to divide global oceans up into seven regions: 
    # Southern Ocean, Indian Ocean, S. Pacific, S. Atlantic, N. Pacific, N. Atlantic, and Arctic

    import iris

    # Define constraints for each region:
    constraint_so = iris.Constraint(latitude=lambda y: -90 < y <= -50)
    constraint_spac = iris.Constraint(latitude=lambda y: -50 < y <= 0, longitude=lambda x: -283.5 <= x <= -60)
    constraint_satl = iris.Constraint(latitude=lambda y: -50 < y <= 0, longitude=lambda x: -60 < x <= 73.5)
    constraint_natl = iris.Constraint(latitude=lambda y: 0 < y <= 60, longitude=lambda x: -60 < x <= 73.5)
    constraint_npac = iris.Constraint(latitude=lambda y: 0 < y <= 60, longitude=lambda x: -283.5 <= x <= -60)
    constraint_arctic = iris.Constraint(latitude=lambda y: 60 < y <= 90)

    var_so = var.extract(constraint_so)
    var_spac = var.extract(constraint_spac)
    var_satl = var.extract(constraint_satl)
    var_npac = var.extract(constraint_npac)
    var_natl = var.extract(constraint_natl)
    var_arctic  = var.extract(constraint_arctic)

    return var_so, var_spac, var_satl, var_npac, var_natl, var_arctic 
