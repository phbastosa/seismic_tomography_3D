#!/bin/bash

wave="../../src/wave/wave.cpp"
utils="../../src/essentials/utils.cpp"
model="../../src/essentials/model.cpp"
geometry="../../src/essentials/geometry.cpp"

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model $geometry $wave acoustic_main.cpp $flags -o acoustic_main.exe

./acoustic_main.exe ../parameters.txt

rm *.o *.exe