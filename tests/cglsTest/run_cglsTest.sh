#!/bin/bash

library="/usr/local/cuda/lib64/"
include="/usr/local/cuda/include/"

cgls="../../src/cgls/cgls.cu"
utils="../../src/essentials/utils.cpp"

pgc++ -acc -fast -ta=tesla,cc60 $cgls $utils cglsTest.cpp -lcublas -lcusparse -lm -I $include -L $library -o cglsTest.exe 

./cglsTest.exe

rm *.exe *.o