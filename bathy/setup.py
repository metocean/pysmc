#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from distutils.core import setup, Extension
import pkgconfig
from subprocess import call

PYTHON_VERSION=3.5


# Config of netcdf point interpolation module
glib_flags=pkgconfig.parse('gmodule-2.0')

bathy_interp_module = Extension('bathy_interp',
                                sources = ["bathy_interp/etopo_interp_centred.c",
                                           "bathy_interp/bathy_interp_module.c"],
                                        libraries = ['netcdf']+glib_flags['libraries']+['m', 'python'+str(PYTHON_VERSION)+'m'],
                                        library_dirs = glib_flags['library_dirs'],
                                        include_dirs = ['/usr/include/python'+str(PYTHON_VERSION)+'m/']+glib_flags['include_dirs'],
                                        language = 'c',
                                        extra_compile_args = ["-fPIC"],
                                        extra_link_args = ['-Wl,--export-dynamic', "-Wl,--no-undefined", "-Wl,--copy-dt-needed-entries"]
)

if __name__ == '__main__':
    setup(name='smc_tools',
          version='0.1.1',
          description='Toolkit for SMC grids',
          author='MetOcean Solutions',
          author_email='ops@metocean.co.nz',
          ext_modules = [bathy_interp_module],
    )

