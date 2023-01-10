#/bin/bash

utils="../../src/essentials/utils.cpp"
model="../../src/essentials/model.cpp"
eikonal="../../src/eikonal/eikonal.cpp"
geometry="../../src/essentials/geometry.cpp"
tomography="../../src/tomography/tomography.cpp"

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model $geometry $eikonal $tomography tomography_main.cpp $flags -o tomography_main.exe

./tomography_main.exe ../parameters.txt

# Clean it up
rm *.o *.exe 
