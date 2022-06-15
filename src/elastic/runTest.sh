#!/bin/bash

model="../essentials/model/model.cpp"
utils="../essentials/utils/utils.cpp"
geom="../essentials/geometry/geometry.cpp"
elastic="elastic.cpp"

flags="-fopenmp -std=c++11 -g -lm -O3"

g++ $inout $utils $model $geom $elastic refractionAnalysis.cpp $flags -o test.exe

./test.exe

rm *.o *.exe
