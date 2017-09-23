from gridgen.generate import create_grid_smc
lonmin, lonmax =  155.0, 180.0
latmin, latmax =  -55.0,-30.0
dlat = dlon = 0.125

create_grid_smc('nz',latmin,latmax,lonmin,lonmax,dlat,dlon)

