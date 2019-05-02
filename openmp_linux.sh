#!/bin/bash

pwd
cd openmp
pwd

rm -rf build
mkdir build
cd build

echo ">>>>>>building the project source code...<<<<<<"

cmake -DCMAKE_BUILD_TYPE=Release ../

sleep 1

cmake -DCMAKE_BUILD_TYPE=Release ../

sleep 1

make -j8

sleep 1

echo
echo
echo ">>>>>>Now use 1 thread as baseline...<<<<<<"
echo
echo

time ./parallel_proj_bin 1 ../data/bunny.off

sleep 1

echo
echo
echo ">>>>>>Now use 8 threads to parallelize...<<<<<<"
echo
echo

time ./parallel_proj_bin 8 ../data/bunny.off