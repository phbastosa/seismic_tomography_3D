#!/bin/bash

model="../../src/essentials/model.cpp"
utils="../../src/essentials/utils.cpp"
geometry="../../src/essentials/geometry.cpp"
acoustic="../../src/acoustic/acoustic.cpp"

flags="-acc -fast -ta=multicore -lm"

pgc++ $utils $model $geometry $acoustic modelingRefractions.cpp $flags -o modelingRefractions.exe

./modelingRefractions.exe

rm *.o *.exe

# python3 refractionAnalysis.py 

# rm outputs/*.bin 
# rm outputs/*.txt
