#!/bin/bash

model="src/essentials/model/model.cpp"
utils="src/essentials/utils/utils.cpp"
geom="src/essentials/geometry/geometry.cpp"
eikonal="src/eikonal/eikonal.cpp"
tomography="src/tomography/tomography.cpp"

main="main/tomography3D_main.cpp"

essentials="$model $utils $geom"
paramFile="inputs/parameters.txt"
flags="-acc -fast -ta=tesla,cc60 -g -std::c++11 -fopenmp -lm"

target="bin/tomography3D.x" 

pgc++ $flags $main $eikonal $tomography $essentials -o $target   

rm *.o

./$target $paramFile

