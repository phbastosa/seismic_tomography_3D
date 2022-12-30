#!/bin/bash

model="../../essentials/model.cpp"
utils="../../essentials/utils.cpp"
geom="../../essentials/geometry.cpp"
eikonal="../../eikonal/eikonal.cpp"

flags="-fopenmp -fast -acc -ta=tesla,cc60 -std=c++11 -g -O3 -lm"

pgc++ $utils $model $geom $eikonal timeTest.cpp $flags -o test.exe

./test.exe

rm *.o *.exe *.bin

