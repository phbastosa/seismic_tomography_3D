#!/bin/bash

flags="-std=c++11 -fast -ta=tesla,cc60 -lm"

pgc++ utils.cpp utilsTest.cpp $flags -o test.exe

./test.exe

rm *.o *.exe