#!/bin/bash

inout="../essentials/inout/inout.cpp"
model="../essentials/model/model.cpp"
utils="../essentials/utils/utils.cpp"
geom="../essentials/geometry/geometry.cpp"
eikonal="eikonal.cpp"

flags="-fopenmp -fast -acc -ta=tesla,cc60 -std=c++11 -g -O3 -lm"

pgc++ $inout $utils $model $geom $eikonal timeTest.cpp $flags -o test.exe

./test.exe

rm *.o *.exe

