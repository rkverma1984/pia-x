#!/bin/bash

## clean start
rm -rf vmd_frames vmd_average
mkdir -p vmd_frames vmd_average

protein="dat"
lsligand=("jjc8091" "jjc8088")
lsmatrix=("angle" "distance" "carbon" "chi1")
iframe_str=1
iframe_end=3


## loops
for iligand in "${lsligand[@]}"
do
for imatrix in "${lsmatrix[@]}"
do

## individual frames
for iframe in `seq ${iframe_str} ${iframe_end}`
do
echo "${protein} & ${iligand}; matirx: ${imatrix}; doing frame ${iframe}"
catdcd -o frame${iframe}.pdb -otype pdb -stype pdb -s inputs/traj/${iligand}.pdb -first ${iframe} -last ${iframe} -dcd inputs/traj/${iligand}.dcd &>/dev/null
vmd -dispdev text -e vmdscripts/${imatrix}.tcl frame${iframe}.pdb frame${iframe}.pdb &>/dev/null
sed 's/\ /,/g' ${imatrix}.dat > vmd_frames/${protein}.${iligand}_${imatrix}_${iligand}_${iframe}.csv
rm -f ${imatrix}.dat frame${iframe}.pdb &>/dev/null
done

## average 
echo "${protein} & ${iligand}; matirx: ${imatrix}; doing average"
catdcd -o frames.dcd -otype dcd -stype pdb -s inputs/traj/${iligand}.pdb -first ${iframe_str} -last ${iframe_end} -dcd inputs/traj/${iligand}.dcd &>/dev/null
vmd -dispdev text -e vmdscripts/${imatrix}.tcl inputs/traj/${iligand}.pdb frames.dcd &>/dev/null
sed 's/\ /,/g' ${imatrix}.dat > vmd_average/${protein}.${iligand}_${imatrix}.csv
rm -f ${imatrix}.dat frames.dcd &>/dev/null

done
done