#!/bin/bash

inout="../essentials/inout/inout.cpp"
model="../essentials/model/model.cpp"
utils="../essentials/utils/utils.cpp"
geom="../essentials/geometry/geometry.cpp"
eik="eikonal.cpp"

flags="-fast -acc -std=c++11 -g  -ta=tesla,cc60 -lm"

pgc++ $inout $utils $model $geom $eik eikonalTest.cpp $flags -o test.exe

./test.exe

rm *.o *.exe

python3 analyticTest.py 

rm pod* fim* *.txt