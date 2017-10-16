import unittest
import scipy.io as sio
from gridgen.pygridgen import pygridgen
import matplotlib.pyplot as plt
import numpy as np


class Test(unittest.TestCase):

    def setUp(self):
        self.lonmin, self.lonmax =  155.0, 180.0
        self.latmin, self.latmax =  -55.0,-30.0
        self.newdlat = self.newdlon = 0.5
        self.matFnm = 'nz.mat'

    def test_bathy(self):
        self.g = pygridgen(self.lonmin, self.lonmax, self.latmin, self.latmax,
                           self.newdlon, self.newdlat)
        self.g.generate_grid(bathy='/source/gridgen/noaa/reference_data/etopo1.nc')

    def test_compare(self):
        vmin=-500*1000
        vmax=0
        # -- load in mat file
        matDict = sio.loadmat(self.matFnm, squeeze_me=True)
        keys = ['dlon', 'dlat', 'lon', 'lat', 'depth', 'm3', 'm4', 'mask_map',
                'sx1', 'sy1']
        for key in keys:
            setattr(self, key, matDict[key])
        del keys, matDict
        old = self.depth
        new = np.loadtxt('test.depth_ascii')/1000
        f, axs = plt.subplots(1, 2, sharey=True, figsize=(12,5))
        m = axs[0].pcolormesh(new, cmap=plt.get_cmap('RdBu_r'), vmin=-1000,
                vmax=0)
        plt.colorbar(m, ax=axs[0])
        plt.title('gridgen')
        m = axs[1].pcolormesh(new, cmap=plt.get_cmap('RdBu_r'), vmin=-1000,
                vmax=0)
        plt.colorbar(m, ax=axs[1])
        plt.title('pygridgen')

        f, axs = plt.subplots(1, 2, sharey=True, figsize=(12,5))
        diff = new - old
        m = axs[0].pcolormesh((diff/abs(old))*100, cmap=plt.get_cmap('RdBu_r'), vmin=-10,
                vmax=10)
        plt.colorbar(m, ax=axs[0])
        plt.title("Percentage Diff")
        m = axs[1].pcolormesh(diff, cmap=plt.get_cmap('RdBu_r'),
                vmin=-2, vmax=2)
        plt.colorbar(m, ax=axs[1])
        plt.title("Absolute Diff")
        plt.show()


if __name__ == "__main__":
    unittest.main()
