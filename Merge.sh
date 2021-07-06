#!/bin/bash
FILESU=(2016*_URaw_variables.nc)
FILESV=(2016*_VRaw_variables.nc)
len=${#FILESU[@]};
echo ${len};

for ((i=0; i<${len}; i++));
do
echo ${FILESU[$i]}
echo ${FILESV[$i]}
cp ${FILESU[$i]} MERGED_${FILESU[$i]}
ncks -A ${FILESV[$i]}  MERGED_${FILESU[$i]}
done
