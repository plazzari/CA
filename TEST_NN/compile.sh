#! /bin/bash

# Configuration dim cDim
Cdim=2

./define_array.py   $Cdim > define_array.h
./allocate_array.py $Cdim > allocate_array.h

ifort  -O0 -g -traceback -fp-stack-check -check bounds -fpe0 GlobalMem.F90 NearNeighbour.f90 main.f90
