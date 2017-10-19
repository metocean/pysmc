import scipy.io as sio
import pandas as pd
import numpy as np
import os

outdir = '../SMCPy/user_polygons'
if not os.path.isdir(outdir):
    os.makedirs(outdir)

matDict = sio.loadmat('/source/gridgen/noaa/reference_data/optional_coastal_polygons.mat',
                      squeeze_me=True)

recs = pd.DataFrame.from_records(matDict['user_bound'])
for item in  recs.iterrows():
    output = np.array([item[1].x, item[1].y]).transpose()
    np.savetxt(os.path.join(outdir, "user_polygon-%s.txt" % item[1].name),
               output)
