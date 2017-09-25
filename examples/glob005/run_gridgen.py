from SMCPy.gridgen import create_grid
lonmin, lonmax =  0.0, 360.0
latmin, latmax =  -75.0,75.0
dlat = dlon = 0.05

create_grid('glb005',latmin,latmax,lonmin,lonmax,dlat,dlon)

