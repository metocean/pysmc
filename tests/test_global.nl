[files]
OutputDir = ./outnz
#NEMOfile = ../bathy/nz-etopo1-0.25.nc
#NEMOfile = ../bathy/glob-etopo1-3.nc
NEMOfile = ../bathy/glob-etopo1-3.nc

[grid]
Rotated = False
gid = glob3-3
urlat=71
lllat=-71
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
