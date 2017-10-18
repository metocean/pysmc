from gridgen.generate import create_grid_smc
lonmin, lonmax =  0.0, 360.0
latmin, latmax =  -75.0,75.0
dlat = dlon = 0.25

create_grid_smc('glb1d',latmin,latmax,lonmin,lonmax,dlat,dlon)

