#!/usr/bin/env python
# -*- coding: utf-8 -*-
# File              : run_gridgen.py
# Author            : Tom Durrant <t.durrant@metocean.co.nz>
# Date              : 03.11.2017
# Last Modified Date: 07.11.2017
# Last Modified By  : Tom Durrant <t.durrant@metocean.co.nz>
from SMCPy.noaa_gridgen import create_grid
lonmin, lonmax =  0.0, 360.0
latmin, latmax =  -75.0,75.0
dlat = dlon = 1.0

create_grid('glb1',latmin,latmax,lonmin,lonmax,dlat,dlon, isglobal=1, overwrite=True)

