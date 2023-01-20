#!/bin/bash

utils="../../src/essentials/utils.cpp"
geometry="../../src/essentials/geometry.cpp"
flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $geometry $utils geometryTest.cpp $flags -o geometryTest.exe

./geometryTest.exe

rm *.o *.exe

python3 verifyGeometryTest.py 

rm outputs/*.txt