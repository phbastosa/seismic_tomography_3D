#!/bin/bash

#########################################################################################
# Eikonal 3D - GISIS/SHELL - Paulo Bastos (iniciado em 10/09/21)
#########################################################################################
# User area 

nx=441; dx=50   #         
ny=441; dy=50   #         
nz=23;  dz=50   #         
                   
vpModel="velocity3D.bin"

z=[1000,1150]
v=[1500,2000]

xsrc=1100
ysrc=1100
zsrc=50

travelTimes="nobleFSM_${nz}x${nx}x${ny}_${dz}x${dx}x${dy}m.bin"
#########################################################################################
# Processing area
#########################################################################################
echo -e "Processing parameters\n"

python3 buildModel.py $nx $ny $nz $dz $v $z $vpModel
echo -e "Model was built..."

g++ eikonal3D.cpp -lm -O3 -o run.exe
./run.exe $nx $ny $nz $dx $dy $dz $xsrc $ysrc $zsrc $vpModel $travelTimes 
rm run.exe

python3 showResults.py $nx $ny $nz $dx $dy $dz $xsrc $ysrc $zsrc $vpModel $travelTimes

