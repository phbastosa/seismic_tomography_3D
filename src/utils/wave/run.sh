#!/bin/bash

flags="-acc -fast -ta=tesla,cc60"

pgc++ ../file_manager/file_manager.cpp acoustic.cpp wave_main.cpp $flags -lm -o wave.exe

./wave.exe

rm *.o