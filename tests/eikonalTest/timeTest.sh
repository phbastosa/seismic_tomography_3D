#!/bin/bash

model="../../src/essentials/model.cpp"
utils="../../src/essentials/utils.cpp"
geom="../../src/essentials/geometry.cpp"
eikonal="../../src/eikonal/eikonal.cpp"

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model $geom $eikonal timeTest.cpp $flags -o test.exe

./test.exe

rm *.o *.exe 

