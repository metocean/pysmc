[files]
OutputDir = ./out025
# NEMOfile =. /source/gridgen/noaa/reference_data/etopo2.nc
#NEMOfile = ../examples/glob05-3/glb05-3.nc
#NEMOfile = ../examples/glb2d/glb2d.nc
#NEMOfile = ../examples/nz/NZ.nc
#NEMOfile = ../bathy/glob-etopo1-0.0625.nc
NEMOfile = ../bathy/glob-etopo1-1.nc

[grid]
gid = split
urlat = -30
lllon = 190
urlon = 320
mergelat=60

[conventions]
xname = lon
yname = lat
zname = z
zscale = 1.0
xyorder = False

[smc]
smctiers = 3
depthlim = 150.0
depthvar = 0.05

[ww3meta]
latlonscale = 1.0
llcrnrscale = 1.0
lsmdepth = -0.1
mindepth = 10.0
blockscale  = 0.1
