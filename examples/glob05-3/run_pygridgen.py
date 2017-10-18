from SMCPy.pybathy import create_grid
lonmin, lonmax =  0.0, 360
latmin, latmax =  -76.0,76.0
dlat = dlon = 0.0625

create_grid('glb025-3',lonmin,lonmax,latmin,latmax,dlat,dlon)

