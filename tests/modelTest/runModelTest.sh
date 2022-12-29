#!/bin/bash

flags="-acc -fast -ta=tesla,cc60 -lm"

pgc++ modelTest.cpp $flags -o modelTest.exe

./modelTest.exe

rm *.exe

python3 verifyModelTest.py

rm outputs/*.bin