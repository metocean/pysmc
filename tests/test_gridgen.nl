[files]
OutputDir = ./outgridgen
# NEMOfile = ./glb1.mat
NEMOfile = ../examples/glb2d/out/glb1d.mat

[grid]
Rotated = False
gid = gridgen
urlat=76
lllat=-76
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
