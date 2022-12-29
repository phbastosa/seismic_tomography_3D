#!/bin/bash

flags="-acc -fast -ta=tesla,cc60 -lm"

pgc++ timeTest.cpp $flags -o test.exe

./test.exe

rm *.o *.exe *.bin

