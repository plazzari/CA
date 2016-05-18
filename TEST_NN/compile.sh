#! /bin/bash

# Configuration dim cDim
Cdim=3

./define_array.py   $Cdim > define_array.h
./allocate_array.py $Cdim > allocate_array.h

ifort  GlobalMem.F90 NearNeighbour.f90 main.f90
