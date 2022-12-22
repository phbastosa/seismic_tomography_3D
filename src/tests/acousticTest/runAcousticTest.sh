#!/bin/bash

model="../../essentials/model.cpp"
utils="../../essentials/utils.cpp"
geometry="../../essentials/geometry.cpp"
acoustic="../../acoustic/acoustic.cpp"

flags=" -ta=multicore -lm"

pgc++ $utils $model $geometry $acoustic acousticTest.cpp $flags -o acousticTest.exe

./acousticTest.exe

rm *.o *.exe

# python3 verifyAcousticTest.py 
