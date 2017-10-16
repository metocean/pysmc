from gridgen.pygridgen import PyGridgen
import cPickle
import numpy as np

class PyBathy(PyGridgen):

    def __init__(self, x1, x2, y1, y2, dx, dy, output='PyBathy.pkl', **kwargs):
        super(PyBathy, self).__init__(x1, x2, y1, y2, dx, dy, **kwargs)
        self.output = output

    def generate_pickle(self):
        self.logger.info("Generating pickle")
        m3 = np.ones_like(self.depth_out)
        m4 = np.ones_like(self.depth_out)
        mask_map = np.ones_like(self.depth_out)
        sx1 = np.zeros_like(self.depth_out)
        sy1 = np.zeros_like(self.depth_out)
        output = {'dlon': self.dx, 'dlat': self.dy, 'lon': self.lons,
                'lat': self.lats, 'depth': self.depth_out.transpose(), 'm3': m3, 'm4': m4,
                'mask_map': mask_map, 'sx1': sx1, 'sy1':sy1}
        cPickle.dump(output, open(self.output, "wb" ))

    def run(self):
        self.generate_grid(bathy='/source/gridgen/noaa/reference_data/etopo1.nc')
        self.generate_pickle()

def create_grid(gid,latmin,latmax,lonmin,lonmax,dlat,dlon):
    pybath = PyBathy(latmin, latmax, lonmin, lonmax, dlat, dlon,
                     output=gid+'.pkl')
    pybath.run()
