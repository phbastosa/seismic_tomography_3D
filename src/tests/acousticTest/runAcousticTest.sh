#!/bin/bash

model="../../essentials/model.cpp"
utils="../../essentials/utils.cpp"
geometry="../../essentials/geometry.cpp"
acoustic="../../acoustic/acoustic.cpp"

flags=" -ta=multicore -lm"

python3 buildModel.py

pgc++ $utils $model $geometry $acoustic acousticTest.cpp $flags -o acousticTest.exe

./acousticTest.exe parametersTest.txt

rm *.o *.exe

python3 verifyAcousticTest.py 
