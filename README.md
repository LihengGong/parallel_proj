# Assignment_1 [![Build Status](https://travis-ci.com/NYUGraphics/assignment1-LihengGong.svg?token=u6yadprNEVSaZeyz72c8&branch=master)](https://travis-ci.com/NYUGraphics/assignment1-LihengGong)
Ray Tracing

mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE=Release ../
make -j8

time ./Assignment1_bin 4 ../data/bunny.off
