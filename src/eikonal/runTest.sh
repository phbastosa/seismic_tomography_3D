#!/bin/bash

inout="../essentials/inout/inout.cpp"
model="../essentials/model/model.cpp"
utils="../essentials/utils/utils.cpp"
geom="../essentials/geometry/geometry.cpp"
eikonal="eikonal.cpp"

flags="-fast -acc -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $inout $utils $model $geom $eikonal eikonalAccuracyStudy.cpp $flags -o test.exe

./test.exe

rm *.o *.exe

# python3 simbgf2022.py 
