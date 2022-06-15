#!/bin/bash

model="../essentials/model/model.cpp"
utils="../essentials/utils/utils.cpp"
geom="../essentials/geometry/geometry.cpp"
elastic="elastic.cpp"

flags="-fopenmp -acc -fast -ta=multicore -std=c++11 -g -lm -O3"

pgc++ $utils $model $geom $elastic elasticTest.cpp $flags -o test.exe

./test.exe

rm *.o *.exe

ximage n1=1001 <seismogram.bin perc=99 &