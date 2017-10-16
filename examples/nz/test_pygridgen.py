import unittest
import scipy.io as sio
from gridgen.pygridgen import pygridgen
import matplotlib.pyplot as plt
import numpy as np
from numpy.testing import assert_array_equal


class Test(unittest.TestCase):

    def setUp(self):
        self.lonmin, self.lonmax =  155.0, 180.0
        self.latmin, self.latmax =  -55.0,-30.0
        self.newdlat = self.newdlon = 0.25
        self.matFnm = 'nz.mat'

    def test_bathy(self):
        self.g = pygridgen(self.lonmin, self.lonmax, self.latmin, self.latmax,
                           self.newdlon, self.newdlat)
        self.g.generate_grid(bathy='/source/gridgen/noaa/reference_data/etopo1.nc')

    def test_compare(self):
        self.load_mat()
        self.load_py()
        self.check_equality()
        self.plot_diff()

    def load_mat(self):
        # -- load in mat file
        self.matDict = sio.loadmat(self.matFnm, squeeze_me=True)
        # keys = ['dlon', 'dlat', 'lon', 'lat', 'depth', 'm3', 'm4', 'mask_map',
                # 'sx1', 'sy1']

    def load_py(self):
        self.pyDict = {}
        self.pyDict['depth'] = np.loadtxt('test.depth_ascii')/1000
        self.pyDict['lat'] = np.loadtxt('test.depth_ascii')/1000
        self.pyDict['lats'] = np.arange(self.latmin, self.latmax + self.newdlat/2., self.newdlat)
        self.pyDict['lons'] = np.arange(self.lonmin, self.lonmax + self.newdlon/2., self.newdlon)

    def plot_diff(self):
        f, axs = plt.subplots(1, 2, sharey=True, figsize=(12,5))
        m = axs[0].pcolormesh(self.matDict['depth'],
                cmap=plt.get_cmap('viridis'), vmin=-4000, vmax=0)
        plt.colorbar(m, ax=axs[0])
        plt.title('gridgen')
        m = axs[1].pcolormesh(self.pyDict['depth'],
                cmap=plt.get_cmap('viridis'), vmin=-4000, vmax=0)
        plt.colorbar(m, ax=axs[1])
        plt.title('pygridgen')

        diff = self.pyDict['depth'] - self.matDict['depth']
        f, axs = plt.subplots(1, 2, sharey=True, figsize=(12,5))
        m = axs[0].pcolormesh((diff/abs(self.matDict['depth']))*100, cmap=plt.get_cmap('RdBu_r'),
                vmin=-10, vmax=10)
        plt.colorbar(m, ax=axs[0])
        plt.title("Percentage Diff")
        m = axs[1].pcolormesh(diff, cmap=plt.get_cmap('RdBu_r'),
                vmin=-2, vmax=2)
        plt.colorbar(m, ax=axs[1])
        plt.title("Absolute Diff")
        plt.show()

    def check_equality(self):
        print "Checking depths"
        assert_array_equal(self.matDict['depth'], self.pyDict['depth'])

    def pickle(self):
        import cPickle
        cPickle.dump(self, open( "test.pkl", "wb" ))

class TestPickle(Test):

    def test_bathy(self):
        super(TestPickle, self).test_bathy()
        self.generate_pickle()

    def generate_pickle(self):
        print "Generating pickle"
        m3 = np.ones_like(self.g.depth_out)
        m4 = np.ones_like(self.g.depth_out)
        mask_map = np.ones_like(self.g.depth_out)
        sx1 = np.zeros_like(self.g.depth_out)
        sy1 = np.zeros_like(self.g.depth_out)
        output = {'dlon': self.g.dx, 'dlat': self.g.dy, 'lon': self.g.lons,
                'lat': self.g.lats, 'depth': self.g.depth_out.transpose(), 'm3': m3, 'm4': m4,
                'mask_map': mask_map, 'sx1': sx1, 'sy1':sy1}
        import cPickle
        cPickle.dump(output, open("test.pkl", "wb" ))

    def load_py(self):
        print "Reading pickle"
        import cPickle
        self.pyDict = cPickle.load(open("test.pkl", "rb"))
        self.pyDict['depth'] = self.pyDict['depth']

if __name__ == "__main__":
    unittest.main()
