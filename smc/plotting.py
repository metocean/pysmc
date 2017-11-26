"""
Functions to create matplotlib PatchCollections of irregularly spaced SMC
data. Uses the mapping package "cartopy" for coordinate projection (see:
http://scitools.org.uk/cartopy/).

If you do not wish yo install Cartopy, then you can comment out the lines
in 'generate_patch_collection' that perform the coordinate transform.

Chris Bunney, UK Met Office
"""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, PolyCollection
import matplotlib.colors as mpl_colors
from collections import Iterable
from cartopy.crs import PlateCarree
import cartopy.crs as ccrs


class SMCPatchCollection(object):
    """ Container class for SMC Patches """
    def __init__(self, pcol, mask, crs):
        self.pcol = pcol        # PatchCollection
        self.mask = mask        # data mask
        self.crs = crs          # coordinate reference system


def load_ww3_txt(fnm):
    ret=[]
    with open(fnm) as f:
        f.readline()
        for line in f.readlines():
            ret += line.split()
    return np.array(ret, dtype='float')

def plot_region(dataset, smc_patch_collection, lon1=None, lat1=None, lon2=None, lat2=None,
        field='hs', tidx=0, clim=[None,None]):

    if not isinstance(dataset, nc.Dataset):
        raise TypeError("Dataset should be a netCDF4 dataset object")

    ## load field:
    print " + Loading field data"
    fld = dataset.variables[field][tidx,:]

    plot_region_fld(fld, smc_patch_collection, lon1=lon1, lat1=lat1, lon2=lon2,
            lat2=lat2, field=field, tidx=tidx, clim=clim)

def plot_region_txt(fnm, smc_patch_collection, lon1=None, lat1=None, lon2=None, lat2=None,
        field='hs', tidx=0, clim=[None,None]):

    ## load field:
    print " + Loading field data"
    fld = load_ww3_txt(fnm)

    plot_region_fld(fld, smc_patch_collection, lon1=lon1, lat1=lat1,
             lon2=lon2, lat2=lat2, field=field, tidx=tidx, clim=clim)

def plot_region_fld(fldall, smc_patch_collection, lon1=None, lat1=None, lon2=None, lat2=None,
        field='hs', tidx=0, clim=[None,None]):

    if not isinstance(smc_patch_collection, SMCPatchCollection):
        raise TypeError("smc_patch_collection should be a SMCPatchCollection object")

    p = smc_patch_collection.pcol
    m = smc_patch_collection.mask
    crs = smc_patch_collection.crs

    fld = fldall[m]

    # get color indices:
    cnorm = mpl_colors.Normalize(clim[0], clim[1])
    cm = plt.cm.get_cmap()
    cols = cm(cnorm(fld))

    # set colors in patch collection:
    p.set_facecolors(cols)
    #p.set_edgecolors('none')
    p.set_edgecolors('gray')

    print " + Adding patches to axes"
    ax = plt.gca()
    ax.add_collection(p)
    ax.coastlines(resolution='50m')

    # -- temporary while testing --#
    plt.gca().set_extent([lon1, lon2, lat1, lat2])
    #plt.gca().set_extent([lon1-10, lon2+10, lat1-10, lat2+10])
    cax = plt.cm.ScalarMappable(cmap=cm)
    cax.set_array(fld)
    plt.colorbar(cax)

def generate_patch_collection_from_cellfille(cellfile, nlev, x0, y0, dx, dy, jshift=0,
    lon1=None, lat1=None, lon2=None, lat2=None, target_crs=ccrs.PlateCarree()):
    """ Generates matplotlib Polygon patches for each SMC grid cell. Stored as
        a PatchCollection for efficient plotting.

        Inputs
        ======

        cellfile, string
            Filename of SMC cell input file to WW3_GRID.

        lon1, lon2, lat1, lat2; integer, optional:
            Corners of area-of-interest. If not provied then whole grid is used.

        target_crs: cartopy.crs instance
            The target coordinate reference system to use. Default to standard
            PlateCarree. Needs to match the coordinate system you plan to use
            to plot the patches.
    """

    with open(cellfile,'r') as fid:
        line = fid.readline()
        hdr = np.fromstring(line, sep=' ', count=5, dtype=np.int)
        ncel = hdr[0]
        print "[INFO] Reading %d cells:" % ncel

        celdat = np.fromfile(fid, sep=' ', dtype=np.int, count=ncel*5).reshape((-1,5))
        print "[INFO] Read cell file"

    ix = celdat[:,0]
    iy = celdat[:,1] + jshift
    cx = celdat[:,2]
    cy = celdat[:,3]

    cfac = 2**(nlev-1)

    dlon1 = dx / cfac
    dlat1 = dy / cfac

    # these are the lat/lons of the cell's SW corner:
    lons = ix * dlon1 # + 0.5*cx*dlon1
    lats = iy * dlat1 # + 0.5*cy*dlat1

    print lons.min(),lons.max();
    print lats.min(),lats.max()
    lons[lons>180] -= 360.0

    print "[INFO] Calculated cell centres"
    print "[INFO] Smallest cell dlon, dlat:",dlon1,dlat1
    #return lons,lats,ix,iy


    # extract region from data:
    if lon1 is not None:
        m1 = np.logical_and(lons > lon1, lons < lon2)
    else:
        m1 = np.ones(lons.shape, dtype=np.bool)

    if lat1 is not None:
        m2 = np.logical_and(lats > lat1, lats < lat2)
    else:
        m2 = np.ones(lons.shape, dtype=np.bool)

    m = np.logical_and(m1, m2)
    del(m1, m2)

    lons = lons[m]
    lats = lats[m]

    dlon = cx[m] * dlon1
    dlat = cy[m] * dlat1

    # calculate four corners of cell in standard (Plate Carree) coorindate reference system:
    x1 = lons
    x2 = lons + dlon
    y1 = lats
    y2 = lats + dlat

    # Transform 4 corners of cells to target coordinate reference system:
    if target_crs is None:
        # don't use cartopy:
        # NOTE: I have not tested this!
        c1 = [x1,y1]
        c2 = [x2,y1]
        c3 = [x2,y2]
        c4 = [x1,y2]
    else:
        # Use cartopy to project the coordinats onto the target reference system
        src_crs = ccrs.PlateCarree()
        c1 = target_crs.transform_points(src_crs,x1,y1)[:,:2]
        c2 = target_crs.transform_points(src_crs,x2,y1)[:,:2]
        c3 = target_crs.transform_points(src_crs,x2,y2)[:,:2]
        c4 = target_crs.transform_points(src_crs,x1,y2)[:,:2]

    # Get the corners into an array of verices of shape (n,4,2). This equates to the
    # x,y point of the 4 corners of each polygon.
    vrts = np.dstack((c1,c2,c3,c4))
    vrts = vrts.swapaxes(2,1)

    # Use vertices to create a PolyCollection:
    p = PolyCollection(vrts)

    return SMCPatchCollection(p, m, target_crs)

def generate_patch_collection(dataset, lon1=None, lat1=None, lon2=None,
        lat2=None, target_crs=ccrs.PlateCarree(), cfacs=None,
        source_crs=None, verbose=0):
    """ Generates matplotlib Polygon patches for each SMC grid cell. Stored as
        a PatchCollection for efficient plotting.

        Inputs
        ======

        dataset, netCDF4 object
            an open netCDF4 dataset object referencing the SMC netCDF file.

        lon1, lon2, lat1, lat2; integer, optional:
            Corners of area-of-interest. If not provied then whole grid is used.

        target_crs: cartopy.crs instance
            The target coordinate reference system to use. Default to standard
            PlateCarree. Needs to match the coordinate system you plan to use
            to plot the patches.
    """


    if not isinstance(dataset, nc.Dataset):
        raise TypeError("Dataset should be a netCDF4 dataset object")


    try:
        lons = dataset.variables['lons'][:]
    except KeyError:
        lons = dataset.variables['longitude'][:]

    lons[lons>180] -= 360 # TODO: Remove this and do longitude alignment
    try:
        lats = dataset.variables['lats'][:]
    except KeyError:
        lats = dataset.variables['latitude'][:]

    cx = dataset.variables['cx'][:]
    cy = dataset.variables['cy'][:]

    if lon1:
        if source_crs:
            lon1,lat1 = source_crs.transform_point(lon1,lat1,ccrs.PlateCarree())
            lon2,lat2 = source_crs.transform_point(lon2,lat2,ccrs.PlateCarree())

    # extract region from data:
    if verbose > 0: print "generate_patch_collection: Extracting region"
    if lon1 is not None:
        m1 = np.logical_and(lons > lon1, lons < lon2)
    else:
        m1 = np.ones(lons.shape, dtype=np.bool)

    if lat1 is not None:
        m2 = np.logical_and(lats > lat1, lats < lat2)
    else:
        m2 = np.ones(lons.shape, dtype=np.bool)

    m = np.logical_and(m1, m2)

    # cfac filter?
    if cfacs != None:
        if verbose > 0: print "Filtering on celfac of %s" % str(cfacs)
        ## Note: filtering only performed on cy
        #m3 = np.logical_or(cx >= max_cfac, cy >= cfac)
        mc = None
        for c in cfacs:
            if verbose > 0: print "Filtering on celfac of %d" % c
            if mc is None:
                mc = cy == c
            else:
                mc = np.logical_or(mc, cy == c)
        m = np.logical_and(m,mc)

    cx = cx[m]
    cy = cy[m]

    lons = lons[m]
    lats = lats[m]
    #dlon = dataset.variables['cx'][:][m] * dataset.base_lon_size
    #dlat = dataset.variables['cy'][:][m] * dataset.base_lat_size
    dlon = cx * dataset.base_lon_size
    dlat = cy * dataset.base_lat_size


    # calculate four corners of cell in standard (Plate Carree) coorindate reference system:
    fac = 0.5 #if isinstance(target_crs, ccrs.PlateCarree) else 0.55
    x1 = lons - fac * dlon
    x2 = lons + fac * dlon
    y1 = lats - fac * dlat
    y2 = lats + fac * dlat

    # Transform 4 corners of cells to target coordinate reference system:
    if target_crs is None:
        # don't use cartopy:
        # NOTE: I have not tested this!
        c1 = [x1,y1]
        c2 = [x2,y1]
        c3 = [x2,y2]
        c4 = [x1,y2]
    else:
        # Use cartopy to project the coordinats onto the target reference system
        if verbose > 0: print "generate_patch_collection: Transforming corner points (%d cells)" % len(x1)
        src_crs = ccrs.PlateCarree() if source_crs is None else source_crs
        c1 = target_crs.transform_points(src_crs,x1,y1)[:,:2]
        c2 = target_crs.transform_points(src_crs,x2,y1)[:,:2]
        c3 = target_crs.transform_points(src_crs,x2,y2)[:,:2]
        c4 = target_crs.transform_points(src_crs,x1,y2)[:,:2]

    # Get the corners into an array of verices of shape (n,4,2). This equates to the
    # x,y point of the 4 corners of each polygon.
    vrts = np.dstack((c1,c2,c3,c4))
    vrts = vrts.swapaxes(2,1)

    # Use vertices to create a PolyCollection:
    if verbose > 0: print "generate_patch_collection: Creating poly collection"
    p = PolyCollection(vrts)

    if verbose > 0: print "generate_patch_collection: Done"
    return SMCPatchCollection(p, m, target_crs)

def generate_color_array(field, cmin=None, cmax=None):
    """ Generate a normalised colour array from a set of data """
    cnorm = mpl_colors.Normalize(cmin, cmax)
    cm = plt.cm.get_cmap()
    cols = cm(cnorm(field))
    return cols

def plot_region_timeseries(dataset, smc_patch_collection, output_fn_template, tidx=None, lon1=None, lat1=None,
        lon2=None, lat2=None, field='hs', clim=[None,None]):

    if not isinstance(dataset, nc.Dataset):
        raise TypeError("Dataset should be a netCDF4 dataset object")

    # get or generate patch collection and input array mask:
    p = smc_patch_collection.pcol
    m = smc_patch_collection.mask
    crs = smc_patch_collection.crs

    # loop over forecast times:
    fld = dataset.variables[field]

    if tidx is None:
        tidx = np.arange(fld.shape[0])

    if not isinstance(tidx, Iterable):
        tidx = [tidx]

    # set up axes:
    plt.figure()
    ax = plt.axes(projection=crs)
    ax.add_collection(p)
    ax.coastlines(resolution='50m')

    p.set_edgecolors('none')
    for t in tidx:
        print " - plotting time idx %d" % t
        cols = generate_color_array(fld[t,:][m])
        p.set_facecolors(cols)

        fn = output_fn_template.replace("<FCHR>","%03d" % t)
        plt.savefig(fn)

    return p, m

def generate_regular_grid_indices(dataset, rlons, rlats):
    """ Use KDSearch tree to generate a set of indices to convert irregular grided SMC
        locations to a regulat lat/lon grid. Nearest neighbour method is used.
        Search radius is set to 4*smallest cell size. """

    from scipy.spatial import cKDTree
    lons = dataset.variables['lons'][:]
    lats = dataset.variables['lats'][:]
    ll = np.vstack([lons,lats]).T
    tree = cKDTree(ll)

    # align lons:
    rlons[rlons < lons.min()] += 360.0
    rlons[rlons > lons.max()] -= 360.0

    # meshgrid the required regular lat/lons:
    rlons,rlats = np.meshgrid(rlons, rlats)

    # use search tree to find nearest neighbours:
    dmin = np.sqrt(dataset.base_lat_size**2 + dataset.base_lat_size**2) * 4
    q = tree.query(np.vstack([rlons.flatten(),rlats.flatten()]).T, distance_upper_bound=dmin)

    # remove entries from query that have a distance of "inf"
    idx = np.where(np.isfinite(q[0]))[0]

    return np.vstack([idx,q[1][idx]]).T

def to_regular_grid(field, indices, gridsize):
    # indices[:,0] = index into gridded array
    # indices[:,1] = index of nearest neighbour in SMC array
    zi = np.ma.masked_all(gridsize)
    ii = np.unravel_index(indices[:,0],gridsize)
    zi[ii] = field[indices[:,1]]

    return zi

def plot_file(fn, **kwargs):
    from ww3tools.core import createMap
    ax = plt.axes(projection=ccrs.PlateCarree())
    d = nc.Dataset(fn, mode='r')
    patches = generate_patch_collection(d,
                verbose=1
            )
    plot_region(d, patches,
                lon1=d.variables['longitude'][:].min(),
                lon2=d.variables['longitude'][:].max(),
                lat1=d.variables['latitude'][:].min(),
                lat2=d.variables['latitude'][:].max(),
                **kwargs)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Plot SMC file")
    parser.add_argument('fn', type=str,
                        help='file path')
    parser.add_argument('-t' '--tidx', type=int, default=0,
                        help='timestep to plot')
    parser.add_argument('-v' '--variable', type=str, default='hs',
                        help='variable to plot')
    args = parser.parse_args()
    plot_file(args.fn, tidx=args.t__tidx, field=args.v__variable)
    plt.show()

if __name__ == "__main__":
    main()
