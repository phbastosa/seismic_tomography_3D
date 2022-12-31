#/bin/bash

utils="../../src/essentials/utils.cpp"
model="../../src/essentials/model.cpp"
geometry="../../src/essentials/geometry.cpp"
eikonal="../../src/eikonal/eikonal.cpp"
tomography="../../src/tomography/tomography.cpp"

flags="-acc -fast -ta=tesla,cc60 -std=c++11 -g -lm"

python3 generateTomoModels.py

pgc++ $utils $model $geometry $eikonal $tomography tomographyTest.cpp $flags -o tomographyTest.exe

./tomographyTest.exe parametersTest.txt

python3 verifyTomographyTest.py

# Clean it up
rm *.o *.exe 
rm outputs/*.bin 
rm outputs/nodesPosition.txt 
rm outputs/shotsPosition.txt 
rm outputs/convergency.txt