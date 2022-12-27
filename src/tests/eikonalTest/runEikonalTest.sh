#!/bin/bash

python3 simbgf2022.py 1

model="../../essentials/model.cpp"
utils="../../essentials/utils.cpp"
geom="../../essentials/geometry.cpp"
eikonal="../../eikonal/eikonal.cpp"

flags="-fopenmp -fast -acc -ta=tesla,cc60 -std=c++11 -g -lm"

pgc++ $utils $model $geom $eikonal eikonalAccuracyStudy.cpp $flags -o eikonalTest.exe

./eikonalTest.exe

rm *.o *.exe

for s in 1 2 3 4
do
    for i in 100 50 25
    do
        if [ $s = "1" ]; then
            mv "outputs/fim_central_${i}m_times_nr5048_shot_${s}.bin" "outputs/fim_5_${i}m.bin"
            mv "outputs/pod_central_${i}m_times_nr5048_shot_${s}.bin" "outputs/pod_5_${i}m.bin"
            mv "outputs/fsm_central_${i}m_times_nr5048_shot_${s}.bin" "outputs/fsm_5_${i}m.bin"
        fi
        
        mv "outputs/fim_extern_${i}m_times_nr5048_shot_${s}.bin" "outputs/fim_${s}_${i}m.bin"    
        mv "outputs/pod_extern_${i}m_times_nr5048_shot_${s}.bin" "outputs/pod_${s}_${i}m.bin"    
        mv "outputs/fsm_extern_${i}m_times_nr5048_shot_${s}.bin" "outputs/fsm_${s}_${i}m.bin"    
    done
done

python3 simbgf2022.py 0

rm outputs/*.bin
rm outputs/*.txt