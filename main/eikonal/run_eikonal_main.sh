#!/bin/bash

model="../../src/essentials/model.cpp"
utils="../../src/essentials/utils.cpp"
geom="../../src/essentials/geometry.cpp"
eikonal="../../src/eikonal/eikonal.cpp"

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $model $utils $geom $eikonal eikonal_main.cpp $flags -o eikonal_main.exe

./eikonal_main.exe ../parameters.txt

rm *.o *.exe