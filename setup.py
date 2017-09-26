#!/usr/bin/env python

from numpy.distutils.core import setup

def ext_configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)
    config.add_extension('SMCPy.fortran.GenCellSide',sources = ['SMCPy/fortran/GenCellSide.pyf','SMCPy/fortran/GenCellSide.f90'])
    return config

k=ext_configuration(top_path='').todict()
k['packages']=['SMCPy', 'SMCPy.fortran']
k['package_dir']={'SMCPy':'./SMCPy'}

setup(name='SMCPy',
      version='0.1',
      description='MSL SMC Grid Creation',
      author='MetOcean Solutions Ltd.',
      author_email='t.durrant@metocean.co.nz',
      url='http://www.metocean.co.nz/',
      **k
      )

