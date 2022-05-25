#!/bin/bash

flags="-std=c++11 -fast -ta=tesla,cc60 -lm"

pgc++ ../inout/inout.cpp model.cpp modelTest.cpp $flags -o test.exe

./test.exe

rm *.o *.exe saltDome3D_expand_test_z128_x378_y378.bin
