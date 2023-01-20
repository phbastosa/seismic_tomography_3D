#/bin/bash

cgls="../../src/cgls/cgls.cu"
utils="../../src/essentials/utils.cpp"
model="../../src/essentials/model.cpp"
eikonal="../../src/eikonal/eikonal.cpp"
geometry="../../src/essentials/geometry.cpp"
tomography="../../src/tomography/tomography.cpp"

library="/usr/local/cuda/lib64/"
include="/usr/local/cuda/include/"

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm -lcusparse -lcublas"

pgc++ $cgls $utils $model $geometry $eikonal $tomography tomography_main.cpp $flags -I $include -L $library -o tomography_main.exe

./tomography_main.exe ../parameters.txt

# Clean it up
rm *.o *.exe 
