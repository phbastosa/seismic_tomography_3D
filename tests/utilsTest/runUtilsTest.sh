#!/bin/bash

flags="-fast -acc -ta=tesla,cc60 -lm"

pgc++ utilsTest.cpp $flags -o utilsTest.exe

./utilsTest.exe

rm *.exe

python3 verifyUtilsTest.py

rm outputs/*.bin