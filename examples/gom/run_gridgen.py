from SMCPy.gridgen import create_grid
lonmin, lonmax =  260, 315
latmin, latmax =  5.0, 32.0
dlat = dlon = 0.025

create_grid('gom',latmin,latmax,lonmin,lonmax,dlat,dlon)

