#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
smc      : python library for plotting and manipulating smc grids
"""

import os
from setuptools import setup, find_packages

import smc
import SMCPy

install_requires = [
    'netCDF4',
    'matplotlib',
    'cartopy',
    ]

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if __name__ == '__main__':
    setup(name = 'pymsl',
          version = smc.__version__,
          description = 'smc',
          author = "Chris Bunney, Tom Durrant",
          author_email = "chris.bunney@metoffice.co",
          maintainer = "Tom Durrant",
          maintainer_email = "t.durrant@metocean.co.nz",
          url = 'https://github.com/metocean/pymsl',
          install_requires=install_requires,
          long_description=read('README.md'),
          packages=['smc'],
         )
    setup(name = 'SMCPy',
          version = SMCPy.__version__,
          description = 'SMCPy',
          author = "Chris Bunney, Tom Durrant",
          author_email = "chris.bunney@metoffice.co",
          maintainer = "Tom Durrant",
          maintainer_email = "t.durrant@metocean.co.nz",
          url = 'https://github.com/metocean/pymsl',
          install_requires=install_requires,
          long_description=read('README.md'),
          packages=['SMCPy'],
         )
