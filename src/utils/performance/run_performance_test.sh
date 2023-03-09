#!/bin/bash

python3 generate_model.py

library="/usr/local/cuda/lib64/"
include="/usr/local/cuda/include/"

regular="../../geometry/regular/regular.cpp"
circular="../../geometry/circular/circular.cpp"

classic="../../eikonal/classic/classic.cpp"
block_FIM="../../eikonal/block_FIM/block_FIM.cpp"
kernel_FIM="../../eikonal/block_FIM/cuda_kernel.cu"
accurate_FSM="../../eikonal/accurate_FSM/accurate_FSM.cpp"

flags="-acc -fast -ta=tesla,cc60 -lm"

trilinear="../interpolation/trilinear.cpp"
file_manager="../file_manager/file_manager.cpp"

pgc++ -fopenmp performance_test.cpp $file_manager $trilinear $regular $circular $classic $block_FIM $kernel_FIM $accurate_FSM $flags -I $include -L $library --diag_suppress bad_macro_redef -o performance.exe

./performance.exe parameters.txt

rm *.o
rm xyz_*
rm *.exe
rm *.bin