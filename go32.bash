#!/bin/bash

if [ `uname` = "Darwin" ]; then
    mexec="/Applications/Mathematica.app/Contents/MacOS/MathKernel "
else
    mexec="/usr/local/bin/math"
fi



for i in 33
do
cp Mathematica/LinearTimeOptimization/LinearTimeOptimization.m tmp/proc
echo " " >> tmp/proc
echo "exploreData[];" >> tmp/proc

chmod 700 tmp/proc
$mexec < tmp/proc > trace$i.txt &
tail -1 tmp/proc
sleep 5;
done
