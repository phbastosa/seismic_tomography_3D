#!/bin/bash

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ geometryTest.cpp $flags -o geometryTest.exe

./geometryTest.exe

rm *.exe

python3 verifyGeometryTest.py 

rm outputs/*.txt