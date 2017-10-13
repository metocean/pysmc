SMC Python Tools
================

Set of python tools to generate SMC grids

Install
--------

Generate python modules from fortran code (required f2py)::

    cd SMCPY/fortran
    ./generate_fpys.sh
    cd ../..

Standard pip install from there, or developer install as follows::

    pip install --user -e .
