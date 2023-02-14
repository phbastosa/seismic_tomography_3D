#!/bin/bash

python3 simbgf2022.py 1

model="../../src/essentials/model.cpp"
utils="../../src/essentials/utils.cpp"
geom="../../src/essentials/geometry.cpp"
eikonal="../../src/eikonal/eikonal.cpp"

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model $geom $eikonal eikonalAccuracyStudy.cpp $flags -o eikonalTest.exe

./eikonalTest.exe

rm *.o *.exe

python3 simbgf2022.py 0
