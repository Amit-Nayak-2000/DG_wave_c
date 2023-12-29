#!/bin/bash

cd outputs
rm -r *.dat
cd ..
rm -r build
mkdir build && cd build
cmake ..
make -j
mpirun -np 8 main