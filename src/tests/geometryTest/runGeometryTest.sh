#!/bin/bash

geometry="../../essentials/geometry.cpp"
flags="-fopenmp -fast -acc -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $geometry geometryTest.cpp $flags -o geometryTest.exe

./geometryTest.exe

rm *.o *.exe

python3 verifyGeometryTest.py 

rm outputs/*.txt