[files]
OutputDir = ./out
# NEMOfile =. /source/gridgen/noaa/reference_data/etopo2.nc
#NEMOfile = ../examples/glob05-3/glb05-3.nc
NEMOfile = ../examples/glb2d/glb2d.nc
#NEMOfile = ../examples/nz/NZ.nc

[grid]
Rotated = False
rlat = 37.5
rlon = 177.5
dx = 1.0
dy = 1.0
lllon = 100
#lllon = 100
urlon = 180
lllat = -50
urlat = -10
#lllon = -8.0
#lllat = -7.0

[conventions]
xname = lon
yname = lat
zname = z
zscale = 1.0
xyorder = False

[smc]
smctiers = 2
depthlim = 150.0
depthvar = 0.05

[ww3meta]
latlonscale = 1.0
llcrnrscale = 1.0
lsmdepth = -0.1
mindepth = 10.0
blockscale  = 0.1
