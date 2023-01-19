#!/bin/bash

utils="../../src/essentials/utils.cpp"
geometry="../../src/essentials/geometry.cpp"

pgc++ -acc -fast -ta=tesla,cc60 $utils $geometry geometry_main.cpp -lm -o geometry_main.exe

./geometry_main.exe ../parameters.txt

rm *.o *.exe
