#!/bin/bash

model="../../essentials/model.cpp"
utils="../../essentials/utils.cpp"
geometry="../../essentials/geometry.cpp"
acoustic="../../acoustic/acoustic.cpp"

flags="-fopenmp -fast -acc -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model $geometry $acoustic acousticTest.cpp $flags -o acousticTest.exe

./acousticTest.exe

rm *.o *.exe

python3 verifyAcousticTest.py 
