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


class NC2SMC(object):

    """Docstring for NC2SMC. """

    def __init__(self, config):
        self.config = config

    def read_config(self):

        print 'Reading namelist...'
        self.myconfig = ConfigParser.RawConfigParser()
        self.myconfig.read(self.config)

        # file locations
        self.workdir = self.myconfig.get("files","OutputDir")
        self.NEMOfile = self.myconfig.get("files","NEMOfile")

        # set input rotation and grid values
        self.rotated = False
        if self.myconfig.has_option("grid","Rotated"):
            if self.myconfig.get("grid","Rotated") == 'True':
                self.rotated = True
                self.rlat = np.float( self.myconfig.get("grid","rlat") )
                self.rlon = np.float( self.myconfig.get("grid","rlon") )
                self.ingrid_lllon = np.float( self.myconfig.get("grid","lllon") )  #longitude sw corner of nemo t-grid
                self.ingrid_lllat = np.float( self.myconfig.get("grid","lllat") )  #latitude sw corner of nemo t-grid
                self.dx = np.float( self.myconfig.get("grid","dx") )
                self.dy = np.float( self.myconfig.get("grid","dy") )

        # grid constraints for limited areas
        if self.myconfig.has_option("grid","llx") and self.myconfig.has_option("grid","lllon"):
            raise Exception("Can only respect llx or lllon, remove one from config")
        if self.myconfig.has_option("grid","llx"):
            self.llx = np.int( self.myconfig.get("grid","llx") )
        else:
            self.llx = 0
        if self.myconfig.has_option("grid","lly") and self.myconfig.has_option("grid","lllat"):
            raise Exception("Can only respect llx or lllon, remove one from config")
        if self.myconfig.has_option("grid","lly"):
            self.lly = np.int( self.myconfig.get("grid","lly") )
        else:
            self.lly = 0
        if self.myconfig.has_option("grid","urx"):
            self.urx = np.int( self.myconfig.get("grid","urx") ) + 1
        else:
            self.urx = None
        if self.myconfig.has_option("grid","ury"):
            self.ury = np.int( self.myconfig.get("grid","ury") ) + 1
        else:
            self.ury = None


        # conventions for reading file
        self.xname = self.myconfig.get("conventions","xname") #name for longitude variable
        self.yname = self.myconfig.get("conventions","yname") #name for latitude variable
        self.zname = self.myconfig.get("conventions","zname") #name for depth variable
        self.zscale = np.float( self.myconfig.get("conventions","zscale") ) #combined scale and pos-neg depth convention
        self.bathyscale = np.abs(self.zscale)
        self.xyorder = self.myconfig.get("conventions","xyorder") #true if bathy variable is x-y array, false if y-x array

# smc information
        if self.myconfig.has_section("smc"):
            smc = True
            self.smctiers = np.int( self.myconfig.get("smc","smctiers") ) #how many tiers in the smc refinement
            self.smcscale = 2.0 ** ( self.smctiers - 1.0 )
            # smc output file names and info for WW3Meta
            self.WW3Cels = 'ww3Cels.dat'
            self.WW3BPs  = 'ww3BPlist.txt'
            self.WW3Meta  = 'ww3meta_SMC.txt'
            self.WW3GDef  = 'ww3.SMC_grid_def'
            self.unitbathy = 30
            self.idlabathy = 3
            self.idfmbathy = 1
            self.smc_dcheck = False
            if self.myconfig.has_option("smc","depthlim"):
               self.smc_dcheck = True
               self.smc_dlim = np.float( self.myconfig.get("smc","depthlim") )
               self.smc_dvar = np.float( self.myconfig.get("smc","depthvar") )

        # ww3 metadata info
        #latlonscale = np.float( self.myconfig.get("ww3meta","latlonscale") ) #
        #scale factor used for lat-lon dx-dy
        #llscale     = np.float( self.myconfig.get("ww3meta","llcrnrscale") ) # scale factor used for lower left corner position
        self.latlonscale = 1.0 # scale factor used for lat-lon dx-dy - best to hardwire this; using 1.0 will remove a lot of confusion in set-up??
        self.llscale     = 1.0 # scale factor used for lat-lon dx-dy - best to hardwire this; using 1.0 will remove a lot of confusion in set-up??
        self.depthlim    = np.float( self.myconfig.get("ww3meta","lsmdepth") ) # depth at which land-sea mask is applied
        self.moddepthmin = np.float( self.myconfig.get("ww3meta","mindepth") ) # minimum model depth
        self.blockscale  = np.float( self.myconfig.get("ww3meta","blockscale") ) # scale factor for blocking information
        self.mindepth_switch = False
        if self.myconfig.has_option("ww3meta","setmindepth"):
            if self.myconfig.get("ww3meta","setmindepth") == 'True':
                self.mindepth_switch = True


    def read_nc(self):
        # read in the NEMO file data and get x-y shape values
        d = nc.Dataset(self.NEMOfile)

        self.depths = d.variables[self.zname]

        self.lat   = d.variables[self.yname]
        self.lon   = d.variables[self.xname]

        if self.myconfig.has_option("grid","lllon"):
            self.llx = np.argmin(np.abs(self.lon[:] - float(self.myconfig.get("grid","lllon"))))
        if self.myconfig.has_option("grid","lllat"):
            self.lly = np.argmin(np.abs(self.lat[:] - float(self.myconfig.get("grid","lllat"))))
        if self.myconfig.has_option("grid","urlon"):
            self.urx = np.argmin(np.abs(self.lon[:] - float(self.myconfig.get("grid","urlon"))))
        if self.myconfig.has_option("grid","urlat"):
            self.ury = np.argmin(np.abs(self.lat[:] - float(self.myconfig.get("grid","urlat"))))

        if self.xyorder == 'True':
            self.ingrid_xpts = np.shape(self.depths[self.llx:self.urx,self.lly:self.ury])[0]
            self.ingrid_ypts = np.shape(self.depths[self.llx:self.urx,self.lly:self.ury])[1]
        else:
            self.ingrid_xpts = np.shape(self.depths[self.lly:self.ury,self.llx:self.urx])[1]
            self.ingrid_ypts = np.shape(self.depths[self.lly:self.ury,self.llx:self.urx])[0]

        # get grid ll corner and dx,dy values from regular grid file
        if not self.rotated:
            if self.xyorder == 'True':
                self.lons, self.lats = np.meshgrid(self.lon, self.lat)
                self.ingrid_lllon = self.lons[self.lly,self.llx]
                self.ingrid_lllat = self.lats[self.lly,self.llx]
                self.dx = lons[self.lly,self.llx+1] - self.ingrid_lllon
                self.dy = lons[self.lly+1,self.llx] - self.ingrid_lllat
            else:
                self.lats, self.lons = np.meshgrid(self.lat, self.lon)
                self.ingrid_lllon = self.lons[self.llx,self.lly]
                self.ingrid_lllat = self.lats[self.llx,self.lly]
                self.dx = self.lons[self.llx+1,self.lly] - self.ingrid_lllon
                self.dy = self.lats[self.llx,self.lly+1] - self.ingrid_lllat


    def extract_region(self):
        # sort out the number of x and y cells to actually use - factor of smcscale
        self.use_xpts = np.int( self.ingrid_xpts - self.ingrid_xpts % self.smcscale )
        self.use_ypts = np.int( self.ingrid_ypts - self.ingrid_ypts % self.smcscale )
        print 'Read in ',self.ingrid_xpts,' xpts; using ',self.use_xpts,' xpts'
        print 'Read in ',self.ingrid_ypts,' ypts; using ',self.use_ypts,' ypts'

        # calculate the output grid sw corner cell centre
        self.lllon = self.ingrid_lllon - self.dx/2.0  + (self.smcscale/2.0) * self.dx
        self.lllat = self.ingrid_lllat - self.dy/2.0  + (self.smcscale/2.0) * self.dy

        # calculate the output grid ne corner cell centre
        self.urlon = self.ingrid_lllon - self.dx/2.0  + (self.use_xpts - self.smcscale/2.0) * self.dx
        self.urlat = self.ingrid_lllat - self.dy/2.0  + (self.use_ypts - self.smcscale/2.0) * self.dy


        print 'Located input grid sw corner cell center at:',self.ingrid_lllon, self.ingrid_lllat,' this gets used for the ww3_grid.inp metadata'
        print 'Located SMC grid sw corner cell center at:',self.lllon,self.lllat,' this gets used for the grid_def metadata'
        print 'Located SMC grid ne corner cell center at:',self.urlon,self.urlat,' this gets used for the boundary point metadata'


        # establish writedepths array using (lat,lon) convention
        if self.xyorder == True:
            self.writedepths = np.rot90(self.depths[self.llx:self.llx+self.use_xpts,self.lly:self.lly+self.use_ypts]) * self.zscale #need to check if this line works properly!
        else:
            self.writedepths = (self.depths[self.lly:self.lly+self.use_ypts,self.llx:self.llx+self.use_xpts] * self.zscale).filled(999)
        print 'Dimensions of grid for analysis: ',np.shape(self.writedepths)

        self.writedepths[self.writedepths>=0] = 999
        if self.mindepth_switch:
            self.writedepths[(self.writedepths>=-1.0*mindepth) & (self.writedepths<0.0)] = -1.0 * mindepth

        self.nx = np.shape(self.writedepths)[1]
        self.ny = np.shape(self.writedepths)[0]

        self.writemask   = np.ones(np.shape(self.writedepths)) * 3.0

        self.writeblock = np.zeros([self.ny*2,self.nx])


    def tier(self):
        """
        routine to tier up the data
        there should be at least 2 tiers in an SMC model
        otherwise its just a regular grid!!!
        """

        smcscli = np.int(self.smcscale)

        # 1st tier establish locations next to land
        print ''
        print 'Analysing Tier 1'
        for lpy in range(self.ny):
            for lpx in range(self.nx):

                if np.mod(lpx,200) == 0 and np.mod(lpy,200) == 0:
                    print 'Analysed %d' %lpx + ' cells in row %d' %lpy

                if self.writedepths[lpy,lpx] == 999:
                    self.writemask[lpy,lpx] = 0
                    if lpx-1 >= 0:
                        if self.writedepths[lpy,lpx-1] != 999:
                            self.writemask[lpy,lpx-1] = 1
                    if lpx+1 < self.nx:
                        if self.writedepths[lpy,lpx+1] != 999:
                            self.writemask[lpy,lpx+1] = 1
                    if lpy-1 >= 0:
                        if self.writedepths[lpy-1,lpx] != 999:
                            self.writemask[lpy-1,lpx] = 1
                    if lpy+1 < self.ny:
                        if self.writedepths[lpy+1,lpx] != 999:
                            self.writemask[lpy+1,lpx] = 1

        # 2nd tier
        # if switched on (smc_dcheck) this will retain type 1 cells for
        # water below a cut-off depth and where cell-cell depth variability
        # is above a threshold
        print ''
        print 'Analysing Tier 2'
        for lpy in range(0,self.ny,2):
            for lpx in range(0,self.nx,2):

                if np.mod(lpx,200) == 0 and np.mod(lpy,200) == 0:
                    print 'Analysed %d' %lpx + ' cells in row %d' %lpy

                if not np.all( self.writemask[lpy:lpy+2,lpx:lpx+2]==0 ):
                    if np.any( self.writemask[lpy:lpy+2,lpx:lpx+2]==1 ):
                       for sly in range(2):
                           for slx in range(2):
                               if self.writemask[lpy+sly,lpx+slx] > 1:
                                   self.writemask[lpy+sly,lpx+slx] = 1
                    elif self.smc_dcheck:
                       # depth based variability
                       if np.any( self.writedepths[lpy:lpy+2,lpx:lpx+2] >= -1*self.smc_dlim ):
                	   dmax = np.max( np.abs(self.writedepths[lpy:lpy+2,lpx:lpx+2]) )
                	   dmin = np.min( np.abs(self.writedepths[lpy:lpy+2,lpx:lpx+2]) )
                	   dmean = np.mean( np.abs(self.writedepths[lpy:lpy+2,lpx:lpx+2]) )
                	   ddep = (dmax - dmin) / dmean
                	   if ddep > self.smc_dvar:
                               self.writemask[lpy:lpy+2,lpx:lpx+2] = 1
                	   else:
                               self.writemask[lpy:lpy+2,lpx:lpx+2] = 2
                       else:
                           self.writemask[lpy:lpy+2,lpx:lpx+2] = 2
                    else:
                       self.writemask[lpy:lpy+2,lpx:lpx+2] = 2

        # third tier
        if 3 <= self.smctiers:
            print ''
            print 'Analysing Tier 3'
            for lpy in range(0,self.ny,4):
        	for lpx in range(0,self.nx,4):

                    if np.mod(lpx,200) == 0 and np.mod(lpy,200) == 0:
                	print 'Analysed %d' %lpx + ' cells in row %d' %lpy

                    # ensure we never go straight from tier 3 to tier 1 by searching over a +/-1 box
                    if lpy-1 >= 0 and lpy+5 < self.ny and lpx-1 >=0 and lpx+5 < self.nx:

                	if np.all( self.writemask[lpy-1:lpy+5,lpx-1:lpx+5]>=2 ):
                	    self.writemask[lpy:lpy+4,lpx:lpx+4] = 3
                	#else:
                	#    print 'rejecting'

        # set border cells at highest tier value
        print ''
        print 'Applying highest tier to border cells'

        self.bplist = []
        lpy = 0
        for lpx in range(smcscli,self.nx,smcscli):
            if np.all( self.writemask[lpy:lpy+(smcscli+1),lpx-1:lpx+(smcscli+1)]>=self.smctiers-1 ):
                self.writemask[lpy:lpy+smcscli,lpx:lpx+smcscli] = self.smctiers
                bcy = self.lllat
                bcx = self.lllon + np.float(lpx) * self.dx
                if [bcx,bcy] not in self.bplist:
                    self.bplist.append( [bcx,bcy] )
        lpy = self.ny
        for lpx in range(smcscli,self.nx,smcscli):
            if np.all( self.writemask[lpy-(smcscli+1):lpy,lpx-1:lpx+(smcscli+1)]>=self.smctiers-1 ):
                self.writemask[lpy-smcscli:lpy,lpx:lpx+smcscli] = self.smctiers
                bcy = self.urlat
                bcx = self.lllon + np.float(lpx) * self.dx
                if [bcx,bcy] not in self.bplist:
                    self.bplist.append( [bcx,bcy] )
        lpx = 0
        for lpy in range(smcscli, self.ny, smcscli):
            if np.all( self.writemask[lpy-1:lpy+(smcscli+1),lpx:lpx+(smcscli+1)]>=self.smctiers-1 ):
                self.writemask[lpy:lpy+smcscli,lpx:lpx+smcscli] = self.smctiers
                bcy = self.lllat + np.float(lpy) * self.dy
                bcx = self.lllon
                if [bcx,bcy] not in self.bplist:
                    self.bplist.append( [bcx,bcy] )
        lpx = self.nx
        for lpy in range(smcscli,self.ny,smcscli):
            if np.all( self.writemask[lpy-1:lpy+(smcscli+1),lpx-(smcscli+1):lpx]>=self.smctiers-1 ):
                self.writemask[lpy:lpy+smcscli,lpx-smcscli:lpx] = self.smctiers
                bcy = self.lllat + np.float(lpy) * self.dy
                bcx = self.urlon
                if [bcx,bcy] not in self.bplist:
                    self.bplist.append( [bcx,bcy] )

        # set corner cells at highest tier value
        if np.all( self.writemask[0:(smcscli+1),0:(smcscli+1)]>=self.smctiers-1 ):
            self.writemask[0:smcscli,0:smcscli] = self.smctiers
            bcy = self.lllat
            bcx = self.lllon
            if [bcx,bcy] not in self.bplist:
                self.bplist.append( [bcx,bcy] )
        if np.all( self.writemask[0:(smcscli+1),-1*(smcscli+1):]>=self.smctiers-1 ):
            self.writemask[0:smcscli,-1*smcscli:] = self.smctiers
            bcy = self.lllat
            bcx = self.urlon
            if [bcx,bcy] not in self.bplist:
                self.bplist.append( [bcx,bcy] )
        if np.all( self.writemask[-1*(smcscli+1):,0:(smcscli+1)]>=self.smctiers-1 ):
            self.writemask[-1*smcscli:,0:smcscli] = self.smctiers
            bcy = self.urlat
            bcx = self.lllon
            if [bcx,bcy] not in self.bplist:
                self.bplist.append( [bcx,bcy] )
        if np.all( self.writemask[-1*(smcscli+1):,-1*(smcscli+1):]>=self.smctiers-1 ):
            self.writemask[-1*smcscli:,-1*smcscli:] = self.smctiers
            bcy = self.urlat
            bcx = self.urlon
            if [bcx,bcy] not in self.bplist:
                self.bplist.append( [bcx,bcy] )


    def create_cells_lists(self):
        # now go through the mask and establish tier lists
        print ''
        print 'Creating cells list for tiers'

        if 3 <= self.smctiers:
            self.tier3=[]
            self.t3x = []
            self.t3y = []

        self.tier2=[]
        self.t2x = []
        self.t2y = []

        self.tier1=[]
        self.t1x = []
        self.t1y = []

        if 3 <= self.smctiers: # three tiers
            for lpy3 in range(0,self.ny,4):
        	for lpx3 in range(0,self.nx,4):

                    if np.all( self.writemask[lpy3:lpy3+4,lpx3:lpx3+4]==3 ) :

                	mydepth = np.mean( self.writedepths[lpy3:lpy3+4,lpx3:lpx3+4] )
                        if mydepth != 'masked':
                            self.tier3.append( [lpx3,lpy3,4,4,mydepth*-1] )
                            self.t3x.append(lpx3)
                            self.t3y.append(lpy3)

                    else:

                	for lp2y in range(0,2):
                            for lp2x in range(0,2):

                        	myx2 = lpx3+2*lp2x
                        	myy2 = lpy3+2*lp2y

                		if np.all( self.writemask[myy2:myy2+2,myx2:myx2+2] == 2 ):

                		    mydepth = np.mean( self.writedepths[myy2:myy2+2,myx2:myx2+2] )
                                    if mydepth != 'masked':
                                        self.tier2.append( [myx2,myy2,2,2,mydepth*-1] )
                                        self.t2x.append(myx2)
                                        self.t2y.append(myy2)

                        	else:

                		    for lpy1 in range(0,2):
                        		for lpx1 in range(0,2):

                        		    myx1 = myx2+lpx1
                        		    myy1 = myy2+lpy1

                			    if self.writemask[myy1,myx1] == 1:

                				mydepth = self.writedepths[myy1,myx1]
                                                if mydepth != 'masked':
                                                    self.tier1.append( [myx1,myy1,1,1,mydepth*-1] )
                                                    self.t1x.append(myx1)
                                                    self.t1y.append(myy1)
        else: # two tiers
            for lpy2 in range(0,self.ny,2):
            	for lpx2 in range(0,self.nx,2):

                    if np.all( self.writemask[lpy2:lpy2+2,lpx2:lpx2+2]==2 ) :

                	mydepth = np.mean( self.writedepths[lpy2:lpy2+2,lpx2:lpx2+2] )
                        if mydepth != 'masked':
                            self.tier2.append( [lpx2,lpy2,2,2,mydepth*-1] )
                            self.t2x.append(lpx2)
                            self.t2y.append(lpy2)

                    else:

                	for lpy1 in range(0,2):
                            for lpx1 in range(0,2):

                        	myx1 = lpx2+lpx1
                        	myy1 = lpy2+lpy1

                		if self.writemask[myy1,myx1] == 1:

                		    mydepth = self.writedepths[myy1,myx1]
                                    if mydepth != 'masked':
                                        self.tier1.append( [myx1,myy1,1,1,mydepth*-1] )
                                        self.t1x.append(myx1)
                                        self.t1y.append(myy1)


    def sort(self):
        # sort lists by y location (2nd value in lists)
        # this is the needed order for processing cell face arrays
        print ''
        print 'Sorting tier lists'
        if 3 <= self.smctiers:
            self.tier3.sort( key=lambda x: int(x[1]) )
        self.tier2.sort( key=lambda x: int(x[1]) )
        self.tier1.sort( key=lambda x: int(x[1]) )

        if 3 <= self.smctiers:
            self.ntc3 = len(self.tier3)
            print 'Number of tier 3 cells: ',self.ntc3
        self.ntc2 = len(self.tier2)
        print 'Number of tier 2 cells: ',self.ntc2
        self.ntc1 = len(self.tier1)
        print 'Number of tier 1 cells: ',self.ntc1
        self.ttotc = self.ntc1 + self.ntc2
        if 3 <= self.smctiers:
            self.ttotc = self.ttotc + self.ntc3


    def plot(self):
        m = plt.pcolormesh(self.writemask)
        if 3 <= self.smctiers:
            plt.scatter(self.t3x,self.t3y,s=1,marker='.',color='k')
        plt.scatter(self.t2x,self.t2y,s=1,marker='.',color='b')
        plt.scatter(self.t1x,self.t1y,s=1,marker='.',color='r')
        plt.colorbar(m)
        plt.show()


    def write_cell(self):
        # write out the cells file
        print ''
        if not os.path.isdir(self.workdir):
            os.makedirs(self.workdir)
        print 'Writing cell info to '+self.workdir+'/'+self.WW3Cels
        with open(self.workdir+'/'+self.WW3Cels,'w') as inp:

            if self.smctiers == 2:
                inp.write( ' %8d %8d %8d' %tuple([self.ttotc,self.ntc1,self.ntc2]) + '  0  0\r\n' )
            elif self.smctiers == 3:
                inp.write( ' %8d %8d %8d %8d' %tuple([self.ttotc,self.ntc1,self.ntc2,self.ntc3]) + '  0\r\n' )

            for lp in range(len(self.tier1)):
                inp.write(' %5d %5d %2d %2d %5d' %tuple(self.tier1[lp]) + '\r\n')

            for lp in range(len(self.tier2)):
                inp.write(' %5d %5d %2d %2d %5d' %tuple(self.tier2[lp]) + '\r\n')

            if 3 <= self.smctiers:
                for lp in range(len(self.tier3)):
                    inp.write(' %5d %5d %2d %2d %5d' %tuple(self.tier3[lp]) + '\r\n')

            inp.close()


    def write_meta(self):
        # write info to metadata file

        # calculating output metadata here
        gdx = self.smcscale * self.dx / self.latlonscale # values for grid_def file - based on largest smc cell size
        gdy = self.smcscale * self.dy / self.latlonscale # values for grid_def file - based on largest smc cell size

        # calculate limits on CFL and 2nd order swell age - based on largest smc cell size
        #maxlat  = (self.lly / self.llscale) + np.float(self.ny) * (self.dy / self.latlonscale)
        maxlat  = max(abs(self.lat[self.lly]), abs(self.lat[self.ury]))
        minlon  = 1853.0 * 60.0 * (self.smcscale * self.dx / self.latlonscale) * np.cos(np.pi*maxlat/180.0)
        maxcg   = 1.4 * 9.81 * 25.0 / (4.0 * np.pi)
        cflstep = minlon / maxcg

        sagemax = 0.5 * minlon**2.0 * 12.0 / ((2.0*np.pi*maxcg/24.0)**2.0 * cflstep)

        # write grid data to grid.inp metadata file
        # note grid parameters are defined by the samllest cell size
        # and use the small cell centre for the sw corner
        print ''
        print 'Writing WW3 metadata to '+self.workdir+'/'+self.WW3Meta
        with open(self.workdir+'/'+self.WW3Meta,'w') as inp:

            inp.write('$ Grid minimum cell dx: %.2f' %minlon +'m at latitude %.3f' %maxlat +' degrees\r\n')
            __import__('ipdb').set_trace()
            inp.write('$ CFL minimum timestep (needs rounding down): %i' %cflstep +' seconds\r\n')
            inp.write('$ Estimated maximum swell age for 24 direction spectrum: %i' %sagemax +' seconds\r\n')
            if self.mindepth_switch:
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
            inp.write(' %i' %self.nx + ' %i' %self.ny +'\r\n')
            inp.write(' %8.6f' %self.dx +' %8.6f' %self.dy +' %5.1f' %self.latlonscale +'\r\n')
            inp.write(' %9.6f' %self.ingrid_lllon +' %9.6f' %self.ingrid_lllat +' %5.1f' %self.llscale +'\r\n')
            inp.write(' %5.2f' %self.depthlim +' %5.1f' %self.moddepthmin +' %i' %self.unitbathy +' %5.1f' %self.bathyscale \
                      +' %i' %self.idlabathy +' %i' %self.idfmbathy +" '(....)' 'NAME' '%s" %self.WW3Cels +"'\r\n")
            inp.write('$\r\n')

            inp.close()


        # write grid data to grid_def file
        # note grid_def parameters for pre-procesing are defined by the largest cell size
        # and use a largest cell centre for the sw corner
        nxdef = np.int( self.nx / self.smcscale )
        nydef = np.int( self.ny / self.smcscale )
        print ''
        print 'Writing grid_def metadata to '+self.workdir+'/'+self.WW3GDef
        with open(self.workdir+'/'+self.WW3GDef,'w') as inp:

            inp.write(' %i' %nxdef + ' %i' %nydef +'\r\n')
            inp.write(' %9.6f' %self.lllon +' %9.6f' %self.lllat +' %8.6f' %gdx +' %8.6f' %gdy +'\r\n')
            # inp.write(' %6.2f' %rlon +' %6.2f' %rlat +'\r\n') # for standard grid set-up

            inp.close()

    def write_bnd(self):
        # write out the boundary points file
        print ''
        print 'Writing boundary point info (real world lat lons) to '+self.workdir+'/'+self.WW3BPs
        with open(self.workdir+'/'+self.WW3BPs,'w') as inpbp:
            for lp in range(len(self.bplist)):
        	writeBP( self.bplist[lp], inpbp, rotated=self.rotated )
            inpbp.close()

    def run(self):
        self.read_config()
        self.read_nc()
        self.extract_region()
        self.tier()
        self.create_cells_lists()
        self.sort()
        self.plot()
        self.write_cell()
        self.write_meta()
        self.write_bnd()



