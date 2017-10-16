from SMCPy.pybathy import create_grid
lonmax, lonmax =  155.0, 180.0
latmin, latmax =  -55.0,-30.0
dlat = dlon = 0.25

create_grid('nz',latmin,latmax,lonmax,lonmax,dlat,dlon)

