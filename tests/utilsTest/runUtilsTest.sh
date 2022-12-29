#!/bin/bash

flags="-fast -acc -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ utilsTest.cpp $flags -o utilsTest.exe

./utilsTest.exe

rm *.exe

python3 verifyUtilsTest.py

rm outputs/*.bin