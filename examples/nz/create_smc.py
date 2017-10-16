# -*- coding: utf-8 -*-
"""
GOMSMCGrid.py

Author: Tom Durrant
Date  : Thu Oct 15 11:59:20 2015
Mod   : Sun Jan 24 2016 [global grid]

Create the SMC2510km grid for Gulf of Mexico (GOM).
"""

# -- m. module
import numpy as np
import scipy.io as sio
import cartopy.crs as ccrs
from datetime import datetime
import matplotlib.pyplot as plt

#import sys; sys.path.append('/media/qliu/Zuoer_dream/mylib')
#import pyutil.SMCGrid as smc
from SMCPy import SMCGrid as smc

# -- 1. important parms
debug = 1
genGrid = 1
matFnm = 'out/nz.mat'
proj=ccrs.Robinson(central_longitude=180.)
#proj=ccrs.Robinson(central_longitude=0.)

# -- 2. gen grid
glbBathy = smc.MatBathy(matFnm, debug=debug)
if genGrid:
    smc.GenSMCGrid(bathy_obj=glbBathy, depth_threshold=150, debug=debug, gen_cell_sides=False)
    # smc.GenSMCGrid(bathy_obj=glbBathy, debug=debug, gen_cell_sides=True)

# -- 3. Create Grid from file and plot the cells
smcFnm = 'NZCell.dat'
obsFnm = 'NZObs.dat'
glbSMC = smc.UnSMC(smcFnm, dlon=glbBathy.dlon, dlat=glbBathy.dlat,
                   refp=(glbBathy.zlon, glbBathy.zlat))#, proj=proj)
glbSMC.readObs(obsFnm)

## -- generate FaceArray and sort [*.d generated by Fortran 90]
#smc.SortFaceArray('GLO1DEGISide.d', 'GOM2510JSide.d')

# -- 3.1 obstruction plot
plot_kws = dict(txtloc=(0.25, 0.88), dotsize=0.01,
                cax_kws=dict(width='2.5%', height='90%', borderpad=.5,
                             loc=6, bbox_to_anchor=(0.975, 0., 1, 1)),
                txtSize=5, cb_kws=dict(orientation='vertical'),
                cbtxtSize=5)

fig, axs = smc.CartopyMap(proj, coast=True, gridbase=30., nrows=2, ncols=1,
                         figsize=(4, 7.5))
glbSMC.genPlot(filled=True, ax=axs[0], plot_var='sx', center=False, **plot_kws)
glbSMC.genPlot(filled=True, ax=axs[1], plot_var='sy', center=False, **plot_kws)
#plt.tight_layout()
plt.savefig(obsFnm[:-4]+'.png')
plt.close()

# -- cell plot
fig, ax = smc.CartopyMap(proj, coast=False, gridbase=45, figsize=(6, 4))
glbSMC.genPlot(filled=False, ax=ax, plot_var='depth', center=True, **plot_kws)

#plt.savefig(smcFnm[:-4]+'.pdf')
plt.savefig(smcFnm[:-4]+'.png')
plt.show()
