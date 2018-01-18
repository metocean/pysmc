import numpy as np
import netCDF4 as nc
import ensemble.smc.regrid as smc_regrid

class SMCIndicies(object):
    """ Storage class for SMC regridding indices """ 

    def __init__(self, nsea):

        ## SMC grid interpolation parameters:
        self.nsea = nsea
        self.xidx = np.ndarray((nsea,), np.int32)
        self.yidx = np.ndarray((nsea,), np.int32)
        self.xspan = np.zeros((nsea,), np.int32)
        self.yspan = np.zeros((nsea,), np.int32)
        self.wts = np.zeros((nsea,), np.float32)

        ## design grid definition:
        self.x0 = None
        self.y0 = None
        self.dx = None
        self.dy = None
        self.nx = None
        self.ny = None

    def __str__(self):
        s = str()
        s+= "Regular grid definition:\n"
        s+= "========================\n"
        s+= "Lons: %.2f -> %.2f (nx=%d, dx=%.4f)\n" % (
                self.x0, self.x0 + (self.nx-1) * self.dx, self.nx, self.dx)
        s+= "Lats: %.2f -> %.2f (ny=%d, dy=%.4f)\n" % (
                self.y0, self.y0 + (self.ny-1) * self.dy, self.ny, self.dy)
        s+= "\n"
        s+= "Interpolation parameters\n"
        s+= "========================\n"
        s+= "Max x/y span: %d/%d\n" % (self.xspan.max(), self.yspan.max())
    
        return s

    def lons(self):
        """ Returns longitudes of cell centres """
        return self.x0 + np.arange(self.nx) * self.dx + 0.5 * self.dx

    def lats(self):
        """ Returns latitudes of cell centres """
        return self.y0 + np.arange(self.ny) * self.dy + 0.5 * self.dy

    def lon_bounds(self):
        """ Returns longitudes boundary of cells (nlon + 1) """
        return self.x0 + np.arange(self.nx+1) * self.dx

    def lat_bounds(self):
        """ Returns latitudes boundary of cells (nlat + 1) """
        return self.y0 + np.arange(self.ny+1) * self.dy

##-- ENDCLASS SMCIndicies


def to_regular_grid_indices_old_format(dataset, sxo, syo, exo, eyo, celfac):
    """ Calculate interpolation indices and weights for regridding the 
        SMC irregular grid back onto a regular lat/lon grid.
    """

    dlon = dataset.base_lon_size
    dlat = dataset.base_lat_size
    x0 = dataset.first_lon
    y0 = dataset.first_lat

    NRLv = 3    # TODO DONT HARDCODE
    cfac = 2**(NRLv - 1)

    sx = dlon * cfac
    sy = dlat * cfac
    print sx,sy

    ilon = dataset.variables['ix'][:]
    ilat = dataset.variables['iy'][:]
    cellx = dataset.variables['cx'][:]
    celly = dataset.variables['cy'][:]


    # SW Corner of grid origin cell:
    cx0 = x0 - sx / 2.
    cy0 = y0 - sy / 2.

    # Get start lat,lon of regular grid (must be aligned with SMC grid edges):
    if np.abs(sxo + 999.9) < 1e-4:
        sxo = cx0 # use SMC grid origin
    else:
        # align to grid cell:
        sxo = cx0 + int((sxo - cx0) / sx) * sx

    if np.abs(syo + 999.9) < 1e-4:
        syo = cy0
    else:
        syo = cy0 + int((syo - cy0) / sy) * sy

    # Last lat/lon (-999.9 denotes use full grid)
    if np.abs(exo + 999.9) < 1e-4:
        exo = cx0 + ilon.max() * dlon
    if np.abs(eyo + 999.9) < 1E-4:
        eyo = cy0 + ilat.max() * dlat

    # Ouput grid cell dx/dy will be integer factor of smallest
    # SMC grid cell size:
    dxo = dlon * celfac
    dyo = dlat * celfac

    # Determine number of cells in grid:
    nxo = int((exo - sxo) / dxo) + 1
    nyo = int((eyo - syo) / dyo) + 1

    print "sxo,syo", sxo,syo
    print "dxo,dyo", dxo,dyo
    print "nxo,nyo", nxo,nyo

    nsea = np.alen(ilon)

    smcind = SMCIndicies(nsea)
    smcind.x0 = sxo
    smcind.y0 = syo
    smcind.dx = dxo
    smcind.dy = dyo
    smcind.nx = nxo
    smcind.ny = nyo
    #xidx = np.ndarray((nsea,), np.int32)
    #yidx = np.ndarray((nsea,), np.int32)
    #xspan = np.zeros((nsea,), np.int32)
    #yspan = np.zeros((nsea,), np.int32)
    #wts = np.zeros((nsea,), np.float32)


    # Loop over cell array and calculate regidding factors:
    for isea in np.arange(nsea):

        # Get grid cell size:
        mx = cellx[isea]
        my = celly[isea]

        # Determine cell lat/lon (this is the SW corner)
        lon = cx0 + ilon[isea] * dlon
        lat = cy0 + ilat[isea] * dlat

        # Align lons
        if lon < sxo:
            lon += 360.
        elif lon > exo:
            lon -= 360.

        # Find first SW cell in design grid:
        ixx = int((lon - sxo) / dxo) # + 1
        iyy = int((lat - syo) / dyo) # + 1

        # range check:
        if ixx < 0 or ixx >= nxo:
            smcind.xidx[isea] = -1
            smcind.yidx[isea] = -1
            continue
       
        if iyy < 0 or iyy >= nyo:
            smcind.xidx[isea] = -1
            smcind.yidx[isea] = -1
            continue


        smcind.xidx[isea] = ixx
        smcind.yidx[isea] = iyy

        # find out how many cells it covers in the x/y directions:
        smcind.xspan[isea] = np.max([1, int(mx / celfac)])
        smcind.yspan[isea] = np.max([1, int(my / celfac)])

        # Do a bit of error checking (non fatal - just produced warning):
        if smcind.xspan[isea] > 1:
            if sxo+ixx*dxo != lon:
                print 'Potential problem with grid cell span:'
                print "Lon,Lat: ",lon,lat
                print "XSPAN,MX:",smcind.xspan[isea], float(mx) / celfac
                print sxo+ixx*dxo,syo+iyy*dyo,dxo,dyo
                return

        # calc cell weight in relation to output grid:
        smcind.wts[isea] = np.min([ 1.0, 
            float(np.min([celfac, mx]) * np.min([celfac, my])) / (celfac**2)
        ])

    ## -- ENDFOR

    return smcind

def to_regular_grid_indices(dataset, sxo, syo, exo, eyo, celfac, NRLv=4):
    """ Calculate interpolation indices and weights for regridding the 
        SMC irregular grid back onto a regular lat/lon grid.

        Works as previous version, but new format files don't have IX/IY
        indices variable - just the lat/lon variable.

    """

    dlon = dataset.base_lon_size
    dlat = dataset.base_lat_size
    x0 = dataset.first_lon
    y0 = dataset.first_lat

    #NRLv = 4    # TODO DONT HARDCODE
    cfac = 2**(NRLv - 1)

    sx = dlon * cfac
    sy = dlat * cfac

    try:
        lons = dataset.variables['longitude'][:]
        lats = dataset.variables['latitude'][:]
    except:
        lons = dataset.variables['lons'][:]
        lats = dataset.variables['lats'][:]

    cellx = dataset.variables['cx'][:]
    celly = dataset.variables['cy'][:]


    # SW Corner of grid origin cell:
    cx0 = x0 - sx / 2.
    cy0 = y0 - sy / 2.

    # Get start lat,lon of regular grid (must be aligned with SMC grid edges):
    if np.abs(sxo + 999.9) < 1e-4:
        sxo = cx0 # use SMC grid origin
    else:
        # align to grid cell:
        sxo = cx0 + int((sxo - cx0) / sx) * sx

    if np.abs(syo + 999.9) < 1e-4:
        syo = cy0
    else:
        syo = cy0 + int((syo - cy0) / sy) * sy

    # Last lat/lon (-999.9 denotes use full grid)
    if np.abs(exo + 999.9) < 1e-4:
        exo = lons.max()
    if np.abs(eyo + 999.9) < 1E-4:
        eyo = lats.max()

    # Ouput grid cell dx/dy will be integer factor of smallest
    # SMC grid cell size:
    dxo = dlon * celfac
    dyo = dlat * celfac

    # Determine number of cells in grid:
    nxo = int((exo - sxo) / dxo) + 1
    nyo = int((eyo - syo) / dyo) + 1

    nsea = np.alen(lons)

    smcind = SMCIndicies(nsea)
    smcind.x0 = sxo
    smcind.y0 = syo
    smcind.dx = dxo
    smcind.dy = dyo
    smcind.nx = nxo
    smcind.ny = nyo
    #xidx = np.ndarray((nsea,), np.int32)
    #yidx = np.ndarray((nsea,), np.int32)
    #xspan = np.zeros((nsea,), np.int32)
    #yspan = np.zeros((nsea,), np.int32)
    #wts = np.zeros((nsea,), np.float32)


    # Loop over cell array and calculate regidding factors:
    for isea in np.arange(nsea):

        # Get grid cell size:
        mx = cellx[isea]
        my = celly[isea]

        # Determine cell lat/lon (this is the SW corner)
        #lon = cx0 + ilon[isea] * dlon
        #lat = cy0 + ilat[isea] * dlat
        lon = lons[isea] - 0.5 * mx * dlon
        lat = lats[isea] - 0.5 * my * dlat

        # Align lons
        if lon < sxo:
            lon += 360.
        elif lon > exo:
            lon -= 360.

        # Find first SW cell in design grid:
        ixx = int((lon - sxo) / dxo) # + 1
        iyy = int((lat - syo) / dyo) # + 1

        # range check:
        if ixx < 0 or ixx >= nxo:
            smcind.xidx[isea] = -1
            smcind.yidx[isea] = -1
            continue
       
        if iyy < 0 or iyy >= nyo:
            smcind.xidx[isea] = -1
            smcind.yidx[isea] = -1
            continue


        smcind.xidx[isea] = ixx
        smcind.yidx[isea] = iyy

        # find out how many cells it covers in the x/y directions:
        smcind.xspan[isea] = np.max([1, int(mx / celfac)])
        smcind.yspan[isea] = np.max([1, int(my / celfac)])
        
    #    print isea,ixx,iyy,smcind.xspan[isea],smcind.yspan[isea]

        # Do a bit of error checking (non fatal - just produced warning):
        if smcind.xspan[isea] > 1:
            if sxo+ixx*dxo != lon:
                print 'Potential problem with grid cell span:'
                print "Lon,Lat: ",lon,lat
                print "XSPAN,MX:",smcind.xspan[isea], float(mx) / celfac
                print sxo+ixx*dxo,syo+iyy*dyo,dxo,dyo
                return

        # calc cell weight in relation to output grid:
        smcind.wts[isea] = np.min([ 1.0, 
            float(np.min([celfac, mx]) * np.min([celfac, my])) / (celfac**2)
        ])

    ## -- ENDFOR

    return smcind

def regular_grid_interp_indicies(dataset, sxo, syo, exo, eyo, celfac, ncel=4):
    """ Interface to C routine """

    _opened_dataset = False
    if not isinstance(dataset, nc.Dataset):
        # assume filenmae passed in, open dataset
        dataset = nc.Dataset(dataset, mode='r')
        _opened_dataset = True
        

    ret = smc_regrid.regrid_indicies(
        dataset.variables['longitude'][:],
        dataset.variables['latitude'][:],
        dataset.variables['cx'][:],
        dataset.variables['cy'][:],
        dataset.first_lon,
        dataset.first_lat,
        dataset.base_lon_size,
        dataset.base_lat_size,
        [sxo, exo, syo, eyo],
        celfac, ncel)

    if _opened_dataset:
        dataset.close()

    smcind = SMCIndicies(ret[0].size)

    smcind.xidx = ret[0]
    smcind.yidx = ret[1]
    smcind.xspan = ret[2]
    smcind.yspan = ret[3]
    smcind.wts = ret[4]
    smcind.x0 = ret[5]
    smcind.y0 = ret[6]
    smcind.dx = ret[7]
    smcind.dy = ret[8]
    smcind.nx = ret[9]
    smcind.ny = ret[10]

    return smcind



def __to_regular_grid_numpy(smcind, dat, fast=True):
    """ Takes a seapoint array of data and regrids onto the target grid
        specified in 'smcind' using precalculatd interpolation indicies.
    """

    undef = -32768.0

    # Initialise coverage and output arrays:
    cov = np.zeros((smcind.nx, smcind.ny), np.float32)
    xy = np.zeros((smcind.nx, smcind.ny), np.float32)

    # get only good cells:
    idx = np.where(smcind.xidx != -1)[0]

    if fast:
        inds = idx
    else:
        inds = np.arange(smcind.nsea)

    for isea in inds:


#        if dat[isea] == undef:
#            continue    # MDI (should not be any in SMC grid)

        
        # Broadcast (is actually slower than looping):
#        ix = np.arange(smcind.xidx[isea], np.min([smcind.xidx[isea] +
#                    smcind.xspan[isea] + 1, smcind.nx]))
#        iy = np.arange(smcind.yidx[isea], np.min([smcind.yidx[isea] +
#                    smcind.yspan[isea] + 1, smcind.ny]))
#        ix,iy = np.meshgrid(ix,iy)
#        xy[ix, iy] = xy[ix, iy] + dat[isea] * smcind.wts[isea]
#        cov[ix, iy] = cov[ix, iy] + smcind.wts[isea]

        # Loop over number of spanned cells:
        for i in np.arange(smcind.xspan[isea]):
            for j in np.arange(smcind.yspan[isea]):
                ix = smcind.xidx[isea] + i
                iy = smcind.yidx[isea] + j

                # Spans outside of grid?
                if ix >= smcind.nx or iy >= smcind.ny:
                    continue

                # Interpolate:
                xy[ix, iy] = xy[ix, iy] + dat[isea] * smcind.wts[isea]

                # Keep track of how much of cell is (wet) covered:
                cov[ix, iy] = cov[ix, iy] + smcind.wts[isea]

    # Create coastline by masking out areas with < 50% coverage:
    xy = np.ma.masked_array(xy, mask=(cov < 0.5))

    # Scale cells on coastlines (those with coverage < 1.0) be normalised
    # against coverage. Without this step, waves at coast will be lower
    # than they should be:
    m = np.logical_and(cov >= 0.5, cov < 1.0)
    xy[m] *=  1.0 / cov[m]

    return xy

def to_regular_grid(dat, smcind, mdi_to_zero=True):
    """ Entry point for fast C regular gridding module

        Inputs:
           dat : a numpy array of shape (nsea) containing the SMC seapoint data
           smcind : SMCIndicies class, returned from 'to_regular_grid_indices'

        Keywords:
            mdi_to_zero: (True/False) If true, any MDI values in the sea point
            array are set to zero, rather than missing data. Useful for
            partitioned data where you dont want large gaps in the data.

        Returns:
            A two-d numpy array of containing data interpolated onto regular
            grid as defined in smcind.
    """

    # Ensure data has MDI value set instead of mask:
    if np.ma.is_masked(dat):
        dat.fill_value = -32767.0 # match expected fill value in smcint routine
        dat = dat.filled()

    regdat = smc_regrid.to_regular_grid(dat, smcind, mdi_to_zero=mdi_to_zero)

    return np.ma.masked_array(regdat, mask=(regdat==-32767.0))


def quick_regrid_plot(dataset, fld, ind, itime=0, proj=None, coast=None):
    """ Regrid the requested smc fld to a regular grid using the supplied
        indicies and plot on a Cartopy axes """

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs

    if proj is None:
        proj = ccrs.PlateCarree()

    dat = dataset.variables[fld][itime,:]
    dat = to_regular_grid(dat, ind)

    ax = plt.axes(projection=proj)
    lons = ind.lon_bounds()
    lats = ind.lat_bounds()
    ax.pcolormesh(lons, lats, dat)
    if coast:
        ax.coastline(resolution=coast)



