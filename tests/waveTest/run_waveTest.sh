#!/bin/bash

wave="../../src/wave/wave.cpp"
model="../../src/essentials/model.cpp"
utils="../../src/essentials/utils.cpp"
geometry="../../src/essentials/geometry.cpp"

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model $geometry $wave waveTest.cpp $flags -o waveTest.exe

./waveTest.exe

rm *.exe *.o

python3 verifyWaveTest.py