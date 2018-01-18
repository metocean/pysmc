"""
Small script to convert the user polygons shipped with gridgen into text files.
These are used to used to mask out given bodies of water from the grid based on
a flag file. A similar functionality is implemented here in SMCPy. The files
produced here are in the repository, this script is just left here in case
these polygons are updated or added to future gridgen releases
"""

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
    np.savetxt(os.path.join(outdir, "user_polygon-%s.txt" % str(int(item[1].name)+1)),
               output)
