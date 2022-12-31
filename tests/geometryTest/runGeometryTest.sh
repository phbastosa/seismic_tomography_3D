#!/bin/bash

geometry="../../src/essentials/geometry.cpp"
flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $geometry geometryTest.cpp $flags -o geometryTest.exe

./geometryTest.exe

rm *.o *.exe

python3 verifyGeometryTest.py 

rm outputs/*.txt