import unittest
import netCDF4 as nc
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from plotting import plot_region, generate_patch_collection


class Test(unittest.TestCase):

    def setUp(self):
        #self.fn = '/source/pysmc/glb05-3_ec_optest2/glob_cell20180102T12.nc'
        self.fn = '/source/pysmc/latest/glob_cell20180103T12.nc'
        self.ds = nc.Dataset(self.fn, mode='r')

    def plot(self):
        self.ax = plt.axes(projection=self.proj)
        patches = generate_patch_collection(self.ds, verbose=1,
                                            lat1=self.lat1,
                                            lat2=self.lat2,
                                            lon1=self.lon1,
                                            lon2=self.lon2,
                                            target_crs=self.proj,
                                            source_crs=ccrs.PlateCarree())
        plot_region(self.ds, patches,
                    lat1=self.lat1,
                    lat2=self.lat2,
                    lon1=self.lon1,
                    lon2=self.lon2,
                    coast_res='110m')

    def test_nz(self):
        self.proj = ccrs.PlateCarree(central_longitude=180)
        self.lat1 = -50
        self.lat2 = -30
        self.lon1 = 160.0
        self.lon2 = 182.0
        self.plot()

    def test_uk(self):
        self.proj = ccrs.PlateCarree()
        self.lat1 = 35
        self.lat2 = 70
        self.lon1 = -15.0
        self.lon2 = 10.0
        self.plot()

    def tearDown(self):
        plt.show()

if __name__ == "__main__":
    unittest.main()
