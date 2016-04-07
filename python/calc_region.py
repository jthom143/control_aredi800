"""
Script to calculate variable for specific regions in ocean
"""

def calc_weddell_sea(var):
    import iris
    
    constraint = iris.Constraint(latitude=lambda y: -85 < y < -60, longitude=lambda x: -60 < x < 0)
    var_ws = var.extract(constraint)

    return var_ws
