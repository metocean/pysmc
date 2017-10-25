SMC Python Tools
================

This repository contains tools to generate and work with Spherical-Multi-Cell
(SMC) grids for WAVEWATCH III.  Code is based in large part from existing bits
and pieces that have been cobbled together into a common library and modified.
Tools provided here are very much a work in progress, and are in various stages
of completeness. 

The main python class::

    SMCPy/         - Main python class for creating SMC grids
    SMCPy/fortran/ - Fortran code for generating face arrays
    SMCPy/matlab/  - Matlab gridgen code for generating bathy
    examples/      - Examples of how this library is used

Other bits and pieces::
    
    smc            - Python code for generating plots from smc flat output files
    idl            - IDL scripts for producing SMC grids
    bathy          - C code to produce bathymetry
    ww3_src        - Alterations to ww3 code 


Install
--------

Generate python modules from fortran code (required f2py)::

    cd SMCPY/fortran
    ./generate_fpys.sh
    cd ../..

Standard pip install from there, or developer install as follows::

    pip install --user -e .

Example
--------

A simple example can be found in examples/nz/,  consisting of a 1 deg base NZ grid refining down to 0.25 based on a refining depth of 250m. 

Base bathymetry was generated from NOAA gridgen code, and is included here so does not need to be produced, however, if you have gridgen installed, it can be replicated by running the script

    python run_gridgen.py

Note that you may need change the paths in SMCPy/matlab/create_grid_smcbase.m to point to the gridgen routines.


The SMC cell file and face arrays can then be produced by running

    python create_smc.py

This will produce the following SMC grid. 

.. image:: examples/nz/NZCell.png

