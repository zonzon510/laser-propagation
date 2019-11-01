#!/bin/bash

# build laser_propagation
cd laser_propagation && \
# mkdir build && \
cd build && \
# cmake .. && \
# make && \

# return to main directory
cd ../../ && \
cd hhgmax
# build the mex file for matlab
CPPFLAGS="-fopenmp -O3 -ansi -std=c++11" LDFLAGS="$CPPFLAGS" mkoctfile -lgomp --mex hhgmax_lewenstein.cpp && \
mv ./hhgmax_lewenstein.mex ./hhgmax_lewenstein.mexa64
# matlab -nodesktop -r example_dipole_response
# ls


