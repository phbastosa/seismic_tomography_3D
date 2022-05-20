#!/bin/bash

flags="-std=c++11 -fast -ta=tesla,cc60 -lm"

pgc++ ../inout/inout.cpp model.cpp modelTest.cpp $flags -o test.exe

./test.exe

rm *.o *.exe
