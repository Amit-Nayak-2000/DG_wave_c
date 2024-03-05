#!/bin/zsh

cd outputs
rm -r *.dat
cd ..
rm -r build
mkdir build && cd build
cmake -DCMAKE_CXX_COMPILER=/usr/local/Cellar/gcc/13.2.0/bin/g++-13 ..
make -j
mpirun -np 4 main