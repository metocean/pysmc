"""
colors.py

Author: Qingxiang Liu
Date  : Nov 18 2015

scripts for colormap generation and application.
"""
# ---------------------------------------------------------
# -- m. modules
import os
import numpy as np
from scipy import interpolate

import matplotlib as mpl
import matplotlib.pyplot as plt

__all__ = ['gen_cmap', 'ct10', 'ct10l', 'ct10m', 'ct20', 'fte']

# ---------------------------------------------------------
# gen_cmap
def gen_cmap(cmap=None, ncols=None, white=None, below=None, above=None,
             reverse=False, beg=None, end=None, extend='neither',
             color=True):
    """
    gen_cmap
    ========

    Create a refined colormap based on specific requirement.
    This function is extracted from cfplot module (Author: Andy Heaps NCAS-CMS)

    Input args:
    -----------
    **cmap** : name of colormap (including colormap from ncl, matplotlib and
        idl guide lib). For example, sequential(`viridis`, `parula`) and diverging
        (`amwg_blueyellowred`)
    **ncols** : # of colors included in the required colormap
        + when extend='neither', ncols = nlevs - 1
        + when extend='min/max', ncols = nlevs
        + when extend='both', ncols = nlevs + 1
    **white** : the rank of the white color in the colormap sequence(an index,
        like a single integer, or integer array / list (It must be 0-based)
    **below** : # of colors below the midpoint of cb
    **above** : # of colors above the midpoint of cb
    **reverse** : whether to reverse the original colorbar, default is False
    **beg** : the begining # of the colormap
    **end** : the end # of the colormap
    **color** : if True return the rgb arrays;  if False return a cmap class
    **extend** : extend mode used in contourf and colorbar; only works if
        color is False
    **

    Output args:
    ------------
    **cmap** : a colormap instance or rgb arrays

    References:
    -----------
    http://climate.ncas.ac.uk/~andy/cfplot_sphinx/_build/html/colour_scales.html
    http://matplotlib.org/examples/color/colormaps_reference.html
    http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

    Examples:
    ---------
    >>> cmap = gen_cmap(cmap='viridis', ncols=8)
    >>> cmap = gen_cmap(cmap='amwg_blueyellowred', ncols=20)
    """

    if hasattr(plt.cm, cmap):
        cmap = plt.get_cmap(cmap)
        rgb = cmap(np.arange(cmap.N))[:, :-1] # get the rgb values [0-1]
    elif cmap == 'avhrr':
        rgb = AVHRR()
    else:
        cmap_path = os.path.join(os.path.dirname(__file__), 'colourmaps',
                                 cmap+'.rgb')

        rgb = np.loadtxt(cmap_path)
        if np.max(rgb) > 1.: rgb = rgb / 255. # [0-1]

    rgb = rgb[beg:end, :]

    # interpolate to a number of colors
    if ncols is not None:
        N = np.shape(rgb)[0]
        x = np.arange(N)
        x_new = (N-1) / float(ncols-1) * np.arange(ncols)
        x_new = np.where(x_new > N-1, N-1, x_new) # important

        #print N, x_new
        f = interpolate.interp1d(x, rgb, axis=0)
        rgb = f(x_new)

    # above and below
    if above is not None or below is not None:
        # mid-point of colour map
        N = np.shape(rgb)[0]
        mid = N / 2

        # below
        lower = mid if below is None else below

        if below == 1:
            x_below = 0
        else:
            x_below = (mid-1) / float(lower-1) * np.arange(lower)

        # above
        upper = mid if above is None else above

        if upper == 1:
            x_above = N - 1
        else:
            x_above = (mid-1) / float(upper-1) * np.arange(upper) + mid

        x_new = np.append(x_below, x_above)
        x_new = np.where(x_new > N-1, N-1, x_new)

        # interpoloate
        x = np.arange(N)
        f = interpolate.interp1d(x, rgb, axis=0)

        rgb = f(x_new)

    # white
    if white is not None:
        rgb[white] = 1.

    # reverse
    if reverse:
        rgb = rgb[::-1]


    # return
    if color:
        return rgb
    else:
       # extend mode
        if extend == 'both':
            cmap = mpl.colors.ListedColormap(rgb[1:-1])
            cmap.set_under(rgb[0])
            cmap.set_over(rgb[-1])
        elif extend == 'min':
            cmap = mpl.colors.ListedColormap(rgb[1:])
            cmap.set_under(rgb[0])
        elif extend == 'max':
            cmap = mpl.colors.ListedColormap(rgb[:-1])
            cmap.set_over(rgb[-1])
        elif extend == 'neither':
            cmap = mpl.colors.ListedColormap(rgb)

        return cmap

#   # levels
#   if levels is not None:
#       norm = mpl.colors.BoundaryNorm(levels, cmap.N)
#       return cmap, norm
#   else:
#       return cmap
#
#    if norm:
#        return mpl.colors.ListedColormap(rgb)
#    else:
#        return rgb

# ---------------------------------------------------------
# beautiful line colors
ct10 = gen_cmap('tableau10')
ct10l = gen_cmap('tableau10_light')
ct10m = gen_cmap('tableau10_medium')
ct20 = gen_cmap('tableau20')

#/-- fivethirtyeight
fte = ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c']
# ---------------------------------------------------------
def AVHRR(m=256):
    """
    AHVRR colormap used by NOAA Coastwatch.
    """

    x = np.arange(0.0, m) / (m - 1)

    xr = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
    rr = [0.5, 1.0, 1.0, 0.5, 0.5, 0.0, 0.5]

    xg = [0.0, 0.4, 0.6, 1.0]
    gg = [0.0, 1.0, 1.0, 0.0]

    xb = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
    bb = [0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 0.5]

    r = np.interp(x, xr, rr)
    g = np.interp(x, xg, gg)
    b = np.interp(x, xb, bb)

    return np.flipud(np.c_[r, g, b])

# @------------------------------------------------------------------- @
