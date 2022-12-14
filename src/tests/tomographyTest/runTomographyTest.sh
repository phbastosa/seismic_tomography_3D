#/bin/bash

utils="../../essentials/utils.cpp"
model="../../essentials/model.cpp"
geometry="../../essentials/geometry.cpp"
eikonal="../../eikonal/eikonal.cpp"
tomography="../../tomography/tomography.cpp"

# flags="-fopenmp -fast -acc -ta=tesla,cc60 -std=c++11 -g -lm"

# python3 generateTomoModels.py

# pgc++ $utils $model $geometry $eikonal $tomography tomographyTest.cpp $flags -o tomographyTest.exe

# ./tomographyTest.exe outputs/parametersTest.txt

# python3 verifyTomographyTest.py

# Clean it up
rm *.o *.exe 
rm outputs/*.bin 
rm outputs/nodesPosition.txt 
rm outputs/shotsPosition.txt 
rm outputs/convergency.txt