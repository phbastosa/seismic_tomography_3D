#!/bin/bash

utils="../src/essentials/utils.cpp"
model="../src/essentials/model.cpp"
geometry="../src/essentials/geometry.cpp"

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model $geometry seismic_modeling.cpp $flags -o seismic_modeling.exe

./seismic_modeling.exe parameters.txt

rm *.o *.exe