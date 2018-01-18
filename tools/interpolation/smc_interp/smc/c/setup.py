from distutils.core import setup, Extension

module1 = Extension('regrid',
                    include_dirs = ['/usr/local/sci/lib/python2.7/site-packages/numpy/core/include/numpy/'],
                    libraries = ['npymath'],
                    library_dirs = ['/usr/local/sci/lib/python2.7/site-packages/numpy/core/lib/'],
                    sources = ['regrid.c'])

setup (name = 'SMCRegrid',
       version = '1.0',
       description = 'This is a demo package',
       author = 'Chris Bunney, UK Met Office',
       author_email = 'chris.bunney@metoffice.gov.uk',
       ext_modules = [module1])
