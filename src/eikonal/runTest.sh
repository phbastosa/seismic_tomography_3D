#!/bin/bash

python3 simbgf2022.py 1

inout="../essentials/inout/inout.cpp"
model="../essentials/model/model.cpp"
utils="../essentials/utils/utils.cpp"
geom="../essentials/geometry/geometry.cpp"
eikonal="eikonal.cpp"

flags="-fast -acc -ta=tesla,cc60 -std=c++11 -g -lm"

# pgc++ $inout $utils $model $geom $eikonal eikonalAccuracyStudy.cpp $flags -o test.exe
g++ -std=c++11 $inout $utils $model $geom $eikonal eikonalAccuracyStudy.cpp -lm -o test.exe

./test.exe

rm *.o *.exe

# for s in 1 2 3 4
# do
#     for i in 100 50 25
#     do
#         if [ $s = "1" ]; then
#             mv "fim_central_${i}m_times_nr5048_shot_${s}.bin" "fim_5_${i}m.bin"
#             mv "pod_central_${i}m_times_nr5048_shot_${s}.bin" "pod_5_${i}m.bin"
#         fi
        
#         mv "fim_extern_${i}m_times_nr5048_shot_${s}.bin" "fim_${s}_${i}m.bin"    
#         mv "pod_extern_${i}m_times_nr5048_shot_${s}.bin" "pod_${s}_${i}m.bin"    
#     done
# done

# python3 simbgf2022.py 0

# rm *.bin *.txt