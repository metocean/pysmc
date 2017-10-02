#!/bin/bash

rm -fr GenCellSide.pyf GenCellSidemodule.c GenCellSide.so
f2py -m GenCellSide -h GenCellSide.pyf GenCellSide.f90 
f2py -c GenCellSide.pyf GenCellSide.f90
