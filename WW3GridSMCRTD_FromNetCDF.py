#!/usr/bin/env python2.7

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import ConfigParser
import xarray
import os

## program to generate SMC grid from NEMO netCDF bathy file

def writeBP( bcs, inpbp, rotated=False ):

    bcx = bcs[0]
    bcy = bcs[1]

    if rotated:
        bplon, bplat = iris.analysis.cartography.unrotate_pole(np.array([bcx]), np.array([bcy]), rlon, rlat)
        bcx = bplon[0]
        bcy = bplat[0]

    if bcx < 0.0:
        bcx = 360.0 + bcx

    #inpbp.write( '%8.3f' %bcx + ',%8.3f' %bcy + '%8.3f' %bcxr + ',%8.3f' %bcyr + ' 0.0 0.0 1\r\n' )
    inpbp.write( '%8.3f' %bcx + ',%8.3f' %bcy + ' 0.0 0.0 1\r\n' )

    return




print 'Reading namelist...'
myconfig = ConfigParser.RawConfigParser()
myconfig.read('./SMCfromNEMO.nl')

# file locations
workdir = myconfig.get("files","OutputDir")
NEMOfile = myconfig.get("files","NEMOfile")

# set input rotation and grid values
rotated = False
if myconfig.has_option("grid","Rotated"):
    if myconfig.get("grid","Rotated") == 'True':
        rotated = True
        rlat = np.float( myconfig.get("grid","rlat") )
        rlon = np.float( myconfig.get("grid","rlon") )
        ingrid_lllon = np.float( myconfig.get("grid","lllon") )  #longitude sw corner of nemo t-grid
        ingrid_lllat = np.float( myconfig.get("grid","lllat") )  #latitude sw corner of nemo t-grid
        dx = np.float( myconfig.get("grid","dx") )
        dy = np.float( myconfig.get("grid","dy") )

# grid constraints for limited areas
if myconfig.has_option("grid","llx"):
    llx = np.int( myconfig.get("grid","llx") )
else:
    llx = 0
if myconfig.has_option("grid","lly"):
    lly = np.int( myconfig.get("grid","lly") )
else:
    lly = 0
if myconfig.has_option("grid","urx"):
    urx = np.int( myconfig.get("grid","urx") ) + 1
else:
    urx = None
if myconfig.has_option("grid","ury"):
    ury = np.int( myconfig.get("grid","ury") ) + 1
else:
    ury = None

# conventions for reading file
xname = myconfig.get("conventions","xname") #name for longitude variable
yname = myconfig.get("conventions","yname") #name for latitude variable
zname = myconfig.get("conventions","zname") #name for depth variable
zscale = np.float( myconfig.get("conventions","zscale") ) #combined scale and pos-neg depth convention
bathyscale = np.abs(zscale)
xyorder = myconfig.get("conventions","xyorder") #true if bathy variable is x-y array, false if y-x array

# smc information
if myconfig.has_section("smc"):
    smc = True
    smctiers = np.int( myconfig.get("smc","smctiers") ) #how many tiers in the smc refinement
    smcscale = 2.0 ** ( smctiers - 1.0 )
    # smc output file names and info for WW3Meta
    WW3Cels = 'ww3Cels.dat'
    WW3BPs  = 'ww3BPlist.txt'
    WW3Meta  = 'ww3meta_SMC.txt'
    WW3GDef  = 'ww3.SMC_grid_def'
    unitbathy = 30
    idlabathy = 3
    idfmbathy = 1
    smc_dcheck = False
    if myconfig.has_option("smc","depthlim"):
       smc_dcheck = True
       smc_dlim = np.float( myconfig.get("smc","depthlim") )
       smc_dvar = np.float( myconfig.get("smc","depthvar") )

# ww3 metadata info
#latlonscale = np.float( myconfig.get("ww3meta","latlonscale") ) # scale factor used for lat-lon dx-dy
#llscale     = np.float( myconfig.get("ww3meta","llcrnrscale") ) # scale factor used for lower left corner position
latlonscale = 1.0 # scale factor used for lat-lon dx-dy - best to hardwire this; using 1.0 will remove a lot of confusion in set-up??
llscale     = 1.0 # scale factor used for lat-lon dx-dy - best to hardwire this; using 1.0 will remove a lot of confusion in set-up??
depthlim    = np.float( myconfig.get("ww3meta","lsmdepth") ) # depth at which land-sea mask is applied
moddepthmin = np.float( myconfig.get("ww3meta","mindepth") ) # minimum model depth
blockscale  = np.float( myconfig.get("ww3meta","blockscale") ) # scale factor for blocking information
mindepth_switch = False
if myconfig.has_option("ww3meta","setmindepth"):
    if myconfig.get("ww3meta","setmindepth") == 'True':
        mindepth_switch = True


#############
# read in the NEMO file data and get x-y shape values
d = nc.Dataset(NEMOfile)

depths = d.variables[zname]
if xyorder == 'True':
    ingrid_xpts = np.shape(depths[llx:urx,lly:ury])[0]
    ingrid_ypts = np.shape(depths[llx:urx,lly:ury])[1]
else:
    ingrid_xpts = np.shape(depths[lly:ury,llx:urx])[1]
    ingrid_ypts = np.shape(depths[lly:ury,llx:urx])[0]

# get grid ll corner and dx,dy values from regular grid file
if not rotated:
    lat   = d.variables[yname]
    lon   = d.variables[xname]
    if xyorder == 'True':
        lons, lats = np.meshgrid(lon, lat)
        ingrid_lllon = lons[lly,llx]
        ingrid_lllat = lats[lly,llx]
        dx = lons[lly,llx+1] - ingrid_lllon
        dy = lons[lly+1,llx] - ingrid_lllat
    else:
        lats, lons = np.meshgrid(lat, lon)
        ingrid_lllon = lons[llx,lly]
        ingrid_lllat = lats[llx,lly]
        dx = lons[llx+1,lly] - ingrid_lllon
        dy = lats[llx,lly+1] - ingrid_lllat


#############
# sort out the number of x and y cells to actually use - factor of smcscale
use_xpts = np.int( ingrid_xpts - ingrid_xpts % smcscale )
use_ypts = np.int( ingrid_ypts - ingrid_ypts % smcscale )
print 'Read in ',ingrid_xpts,' xpts; using ',use_xpts,' xpts'
print 'Read in ',ingrid_ypts,' ypts; using ',use_ypts,' ypts'

# calculate the output grid sw corner cell centre
lllon = ingrid_lllon - dx/2.0  + (smcscale/2.0) * dx
lllat = ingrid_lllat - dy/2.0  + (smcscale/2.0) * dy

# calculate the output grid ne corner cell centre
urlon = ingrid_lllon - dx/2.0  + (use_xpts - smcscale/2.0) * dx
urlat = ingrid_lllat - dy/2.0  + (use_ypts - smcscale/2.0) * dy


print 'Located input grid sw corner cell center at:',ingrid_lllon,ingrid_lllat,' this gets used for the ww3_grid.inp metadata'
print 'Located SMC grid sw corner cell center at:',lllon,lllat,' this gets used for the grid_def metadata'
print 'Located SMC grid ne corner cell center at:',urlon,urlat,' this gets used for the boundary point metadata'


##############
# prepare the NEMO data for SMC processing
# establish writedepths array using (lat,lon) convention
if xyorder == True:
    writedepths = np.rot90(depths[llx:llx+use_xpts,lly:lly+use_ypts]) * zscale #need to check if this line works properly!
else:
    writedepths = depths[lly:lly+use_ypts,llx:llx+use_xpts] * zscale
print 'Dimensions of grid for analysis: ',np.shape(writedepths)

writedepths[writedepths>=0] = 999
if mindepth_switch:
    writedepths[(writedepths>=-1.0*mindepth) & (writedepths<0.0)] = -1.0 * mindepth

nx = np.shape(writedepths)[1]
ny = np.shape(writedepths)[0]

writemask   = np.ones(np.shape(writedepths)) * 3.0

writeblock = np.zeros([ny*2,nx])


###############
# routine to tier up the data
# there should be at least 2 tiers in an SMC model
# otherwise its just a regular grid!!!

smcscli = np.int(smcscale)

# 1st tier establish locations next to land
print ''
print 'Analysing Tier 1'
for lpy in range(ny):
    for lpx in range(nx):

        if np.mod(lpx,200) == 0 and np.mod(lpy,200) == 0:
            print 'Analysed %d' %lpx + ' cells in row %d' %lpy

        if writedepths[lpy,lpx] == 999:
            writemask[lpy,lpx] = 0
            if lpx-1 >= 0:
                if writedepths[lpy,lpx-1] != 999:
                    writemask[lpy,lpx-1] = 1
            if lpx+1 < nx:
                if writedepths[lpy,lpx+1] != 999:
                    writemask[lpy,lpx+1] = 1
            if lpy-1 >= 0:
                if writedepths[lpy-1,lpx] != 999:
                    writemask[lpy-1,lpx] = 1
            if lpy+1 < ny:
                if writedepths[lpy+1,lpx] != 999:
                    writemask[lpy+1,lpx] = 1

# 2nd tier
# if switched on (smc_dcheck) this will retain type 1 cells for
# water below a cut-off depth and where cell-cell depth variability
# is above a threshold
print ''
print 'Analysing Tier 2'
for lpy in range(0,ny,2):
    for lpx in range(0,nx,2):

        if np.mod(lpx,200) == 0 and np.mod(lpy,200) == 0:
            print 'Analysed %d' %lpx + ' cells in row %d' %lpy

        if not np.all( writemask[lpy:lpy+2,lpx:lpx+2]==0 ):
            if np.any( writemask[lpy:lpy+2,lpx:lpx+2]==1 ):
               for sly in range(2):
                   for slx in range(2):
                       if writemask[lpy+sly,lpx+slx] > 1:
                           writemask[lpy+sly,lpx+slx] = 1
            elif smc_dcheck:
               # depth based variability
               if np.any( writedepths[lpy:lpy+2,lpx:lpx+2] >= -1*smc_dlim ):
        	   dmax = np.max( np.abs(writedepths[lpy:lpy+2,lpx:lpx+2]) )
        	   dmin = np.min( np.abs(writedepths[lpy:lpy+2,lpx:lpx+2]) )
        	   dmean = np.mean( np.abs(writedepths[lpy:lpy+2,lpx:lpx+2]) )
        	   ddep = (dmax - dmin) / dmean
        	   if ddep > smc_dvar:
                       writemask[lpy:lpy+2,lpx:lpx+2] = 1
        	   else:
                       writemask[lpy:lpy+2,lpx:lpx+2] = 2
               else:
                   writemask[lpy:lpy+2,lpx:lpx+2] = 2
            else:
               writemask[lpy:lpy+2,lpx:lpx+2] = 2

# third tier
if 3 <= smctiers:
    print ''
    print 'Analysing Tier 3'
    for lpy in range(0,ny,4):
	for lpx in range(0,nx,4):

            if np.mod(lpx,200) == 0 and np.mod(lpy,200) == 0:
        	print 'Analysed %d' %lpx + ' cells in row %d' %lpy

            # ensure we never go straight from tier 3 to tier 1 by searching over a +/-1 box
            if lpy-1 >= 0 and lpy+5 < ny and lpx-1 >=0 and lpx+5 < nx:

        	if np.all( writemask[lpy-1:lpy+5,lpx-1:lpx+5]>=2 ):
        	    writemask[lpy:lpy+4,lpx:lpx+4] = 3
        	#else:
        	#    print 'rejecting'

# set border cells at highest tier value
print ''
print 'Applying highest tier to border cells'

bplist = []
lpy = 0
for lpx in range(smcscli,nx,smcscli):
    if np.all( writemask[lpy:lpy+(smcscli+1),lpx-1:lpx+(smcscli+1)]>=smctiers-1 ):
        writemask[lpy:lpy+smcscli,lpx:lpx+smcscli] = smctiers
        bcy = lllat
        bcx = lllon + np.float(lpx) * dx
        if [bcx,bcy] not in bplist:
            bplist.append( [bcx,bcy] )
lpy = ny
for lpx in range(smcscli,nx,smcscli):
    if np.all( writemask[lpy-(smcscli+1):lpy,lpx-1:lpx+(smcscli+1)]>=smctiers-1 ):
        writemask[lpy-smcscli:lpy,lpx:lpx+smcscli] = smctiers
        bcy = urlat
        bcx = lllon + np.float(lpx) * dx
        if [bcx,bcy] not in bplist:
            bplist.append( [bcx,bcy] )
lpx = 0
for lpy in range(smcscli, ny, smcscli):
    if np.all( writemask[lpy-1:lpy+(smcscli+1),lpx:lpx+(smcscli+1)]>=smctiers-1 ):
        writemask[lpy:lpy+smcscli,lpx:lpx+smcscli] = smctiers
        bcy = lllat + np.float(lpy) * dy
        bcx = lllon
        if [bcx,bcy] not in bplist:
            bplist.append( [bcx,bcy] )
lpx = nx
for lpy in range(smcscli,ny,smcscli):
    if np.all( writemask[lpy-1:lpy+(smcscli+1),lpx-(smcscli+1):lpx]>=smctiers-1 ):
        writemask[lpy:lpy+smcscli,lpx-smcscli:lpx] = smctiers
        bcy = lllat + np.float(lpy) * dy
        bcx = urlon
        if [bcx,bcy] not in bplist:
            bplist.append( [bcx,bcy] )

# set corner cells at highest tier value
if np.all( writemask[0:(smcscli+1),0:(smcscli+1)]>=smctiers-1 ):
    writemask[0:smcscli,0:smcscli] = smctiers
    bcy = lllat
    bcx = lllon
    if [bcx,bcy] not in bplist:
        bplist.append( [bcx,bcy] )
if np.all( writemask[0:(smcscli+1),-1*(smcscli+1):]>=smctiers-1 ):
    writemask[0:smcscli,-1*smcscli:] = smctiers
    bcy = lllat
    bcx = urlon
    if [bcx,bcy] not in bplist:
        bplist.append( [bcx,bcy] )
if np.all( writemask[-1*(smcscli+1):,0:(smcscli+1)]>=smctiers-1 ):
    writemask[-1*smcscli:,0:smcscli] = smctiers
    bcy = urlat
    bcx = lllon
    if [bcx,bcy] not in bplist:
        bplist.append( [bcx,bcy] )
if np.all( writemask[-1*(smcscli+1):,-1*(smcscli+1):]>=smctiers-1 ):
    writemask[-1*smcscli:,-1*smcscli:] = smctiers
    bcy = urlat
    bcx = urlon
    if [bcx,bcy] not in bplist:
        bplist.append( [bcx,bcy] )


# now go through the mask and establish tier lists
print ''
print 'Creating cells list for tiers'

if 3 <= smctiers:
    tier3=[]
    t3x = []
    t3y = []

tier2=[]
t2x = []
t2y = []

tier1=[]
t1x = []
t1y = []

if 3 <= smctiers: # three tiers
    for lpy3 in range(0,ny,4):
	for lpx3 in range(0,nx,4):

            if np.all( writemask[lpy3:lpy3+4,lpx3:lpx3+4]==3 ) :

        	mydepth = np.mean( writedepths[lpy3:lpy3+4,lpx3:lpx3+4] )
                if mydepth != 'masked':
                    tier3.append( [lpx3,lpy3,4,4,mydepth*-1] )
                    t3x.append(lpx3)
                    t3y.append(lpy3)

            else:

        	for lp2y in range(0,2):
                    for lp2x in range(0,2):

                	myx2 = lpx3+2*lp2x
                	myy2 = lpy3+2*lp2y

        		if np.all( writemask[myy2:myy2+2,myx2:myx2+2] == 2 ):

        		    mydepth = np.mean( writedepths[myy2:myy2+2,myx2:myx2+2] )
                            if mydepth != 'masked':
                                tier2.append( [myx2,myy2,2,2,mydepth*-1] )
                                t2x.append(myx2)
                                t2y.append(myy2)

                	else:

        		    for lpy1 in range(0,2):
                		for lpx1 in range(0,2):

                		    myx1 = myx2+lpx1
                		    myy1 = myy2+lpy1

        			    if writemask[myy1,myx1] == 1:

        				mydepth = writedepths[myy1,myx1]
                                        if mydepth != 'masked':
                                            tier1.append( [myx1,myy1,1,1,mydepth*-1] )
                                            t1x.append(myx1)
                                            t1y.append(myy1)
else: # two tiers
    for lpy2 in range(0,ny,2):
    	for lpx2 in range(0,nx,2):

            if np.all( writemask[lpy2:lpy2+2,lpx2:lpx2+2]==2 ) :

        	mydepth = np.mean( writedepths[lpy2:lpy2+2,lpx2:lpx2+2] )
                if mydepth != 'masked':
                    tier2.append( [lpx2,lpy2,2,2,mydepth*-1] )
                    t2x.append(lpx2)
                    t2y.append(lpy2)

            else:

        	for lpy1 in range(0,2):
                    for lpx1 in range(0,2):

                	myx1 = lpx2+lpx1
                	myy1 = lpy2+lpy1

        		if writemask[myy1,myx1] == 1:

        		    mydepth = writedepths[myy1,myx1]
                            if mydepth != 'masked':
                                tier1.append( [myx1,myy1,1,1,mydepth*-1] )
                                t1x.append(myx1)
                                t1y.append(myy1)


# sort lists by y location (2nd value in lists)
# this is the needed order for processing cell face arrays
print ''
print 'Sorting tier lists'
if 3 <= smctiers:
    tier3.sort( key=lambda x: int(x[1]) )
tier2.sort( key=lambda x: int(x[1]) )
tier1.sort( key=lambda x: int(x[1]) )


if 3 <= smctiers:
    ntc3 = len(tier3)
    print 'Number of tier 3 cells: ',ntc3
ntc2 = len(tier2)
print 'Number of tier 2 cells: ',ntc2
ntc1 = len(tier1)
print 'Number of tier 1 cells: ',ntc1
ttotc = ntc1 + ntc2
if 3 <= smctiers:
    ttotc = ttotc + ntc3


########################
# plot the data in order to check the processing

plt.pcolormesh(writemask)
if 3 <= smctiers:
    plt.scatter(t3x,t3y,s=1,marker='.',color='k')
plt.scatter(t2x,t2y,s=1,marker='.',color='b')
plt.scatter(t1x,t1y,s=1,marker='.',color='r')
#plt.colorbar()
plt.show()


########################
# write out the data

# write out the cells file
print ''
if not os.path.isdir(workdir):
    os.makedirs(workdir)
print 'Writing cell info to '+workdir+'/'+WW3Cels
with open(workdir+'/'+WW3Cels,'w') as inp:

    if smctiers == 2:
        inp.write( ' %8d %8d %8d' %tuple([ttotc,ntc1,ntc2]) + '  0  0\r\n' )
    elif smctiers == 3:
        inp.write( ' %8d %8d %8d %8d' %tuple([ttotc,ntc1,ntc2,ntc3]) + '  0\r\n' )

    for lp in range(len(tier1)):
        inp.write(' %5d %5d %2d %2d %5d' %tuple(tier1[lp]) + '\r\n')

    for lp in range(len(tier2)):
        inp.write(' %5d %5d %2d %2d %5d' %tuple(tier2[lp]) + '\r\n')

    if 3 <= smctiers:
	for lp in range(len(tier3)):
            inp.write(' %5d %5d %2d %2d %5d' %tuple(tier3[lp]) + '\r\n')

    inp.close()


# write info to metadata file

# calculating output metadata here
gdx = smcscale * dx / latlonscale # values for grid_def file - based on largest smc cell size
gdy = smcscale * dy / latlonscale # values for grid_def file - based on largest smc cell size

# calculate limits on CFL and 2nd order swell age - based on largest smc cell size
maxlat  = (lly / llscale) + np.float(ny) * (dy / latlonscale)
minlon  = 1853.0 * 60.0 * (smcscale * dx / latlonscale) * np.cos(np.pi*maxlat/180.0)
maxcg   = 1.4 * 9.81 * 25.0 / (4.0 * np.pi)
cflstep = minlon / maxcg

sagemax = 0.5 * minlon**2.0 * 12.0 / ((2.0*np.pi*maxcg/24.0)**2.0 * cflstep)

# write grid data to grid.inp metadata file
# note grid parameters are defined by the samllest cell size
# and use the small cell centre for the sw corner
print ''
print 'Writing WW3 metadata to '+workdir+'/'+WW3Meta
with open(workdir+'/'+WW3Meta,'w') as inp:

    inp.write('$ Grid minimum cell dx: %.2f' %minlon +'m at latitude %.3f' %maxlat +' degrees\r\n')
    inp.write('$ CFL minimum timestep (needs rounding down): %i' %cflstep +' seconds\r\n')
    inp.write('$ Estimated maximum swell age for 24 direction spectrum: %i' %sagemax +' seconds\r\n')
    if mindepth_switch:
        inp.write('$ Minimum depth set for model at %f' %mindepth + 'm\r\n')
    inp.write('$\r\n')
    inp.write('$ Define grid rules -------------------------------------------------- $\r\n')
    inp.write('$ Four records containing :\r\n')
    inp.write('$  1 NX, NY. As the outer grid lines are always defined as land\r\n')
    inp.write('$    points, the minimum size is 3x3.\r\n')
    inp.write('$  2 Grid increments SX, SY (degr.or m) and scaling (division) factor.\r\n')
    inp.write('$    If NX*SX = 360., latitudinal closure is applied.\r\n')
    inp.write('$  3 Coordinates of (1,1) (degr.) and scaling (division) factor.\r\n')
    inp.write('$  4 Limiting bottom depth (m) to discriminate between land and sea\r\n')
    inp.write('$    points, minimum water depth (m) as allowed in model, unit number\r\n')
    inp.write('$    of file with bottom depths, scale factor for bottom depths (mult.),\r\n')
    inp.write('$    IDLA, IDFM, format for formatted read, FROM and filename.\r\n')
    inp.write('$\r\n')
    inp.write('$ Define grid -------------------------------------------------------- $\r\n')
    inp.write('$\r\n')
    inp.write(' %i' %nx + ' %i' %ny +'\r\n')
    inp.write(' %8.6f' %dx +' %8.6f' %dy +' %5.1f' %latlonscale +'\r\n')
    inp.write(' %9.6f' %ingrid_lllon +' %9.6f' %ingrid_lllat +' %5.1f' %llscale +'\r\n')
    inp.write(' %5.2f' %depthlim +' %5.1f' %moddepthmin +' %i' %unitbathy +' %5.1f' %bathyscale \
              +' %i' %idlabathy +' %i' %idfmbathy +" '(....)' 'NAME' '%s" %WW3Cels +"'\r\n")
    inp.write('$\r\n')

    inp.close()


# write grid data to grid_def file
# note grid_def parameters for pre-procesing are defined by the largest cell size
# and use a largest cell centre for the sw corner
nxdef = np.int( nx / smcscale )
nydef = np.int( ny / smcscale )
print ''
print 'Writing grid_def metadata to '+workdir+'/'+WW3GDef
with open(workdir+'/'+WW3GDef,'w') as inp:

    inp.write(' %i' %nxdef + ' %i' %nydef +'\r\n')
    inp.write(' %9.6f' %lllon +' %9.6f' %lllat +' %8.6f' %gdx +' %8.6f' %gdy +'\r\n')
    # inp.write(' %6.2f' %rlon +' %6.2f' %rlat +'\r\n') # for standard grid set-up

    inp.close()

# write out the boundary points file
print ''
print 'Writing boundary point info (real world lat lons) to '+workdir+'/'+WW3BPs
with open(workdir+'/'+WW3BPs,'w') as inpbp:
    for lp in range(len(bplist)):
	writeBP( bplist[lp], inpbp, rotated=rotated )
    inpbp.close()

