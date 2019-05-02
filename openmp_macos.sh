#!/bin/bash

pwd
cd OpenMP
pwd

rm -rf build
mkdir build
cd build

echo ">>>>>>building the project source code...<<<<<<"

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="/usr/local/opt/llvm/bin/clang" -DCMAKE_CXX_COMPILER="/usr/local/opt/llvm/bin/clang++" ../

sleep 1

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="/usr/local/opt/llvm/bin/clang" -DCMAKE_CXX_COMPILER="/usr/local/opt/llvm/bin/clang++" ../

sleep 1

make -j8

sleep 1

echo
echo
echo ">>>>>>First run in single thread as the base line...<<<<<<"
echo
echo

time ./parallel_proj_bin 1 ../data/bunny.off

echo
echo
echo ">>>>>>Now run in 8 threads...<<<<<<"
echo
echo

time ./parallel_proj_bin 8 ../data/bunny.off