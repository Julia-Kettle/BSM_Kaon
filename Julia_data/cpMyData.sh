#!/bin/bash

for lat in 48cubed 48cubedfine 64cubed
do
for c in 0 2 4 6 8
do
for f in /Home/s1035546/Fits/${lat}/bag/channel${c}/data/*.dat
do

name=`basename $f`
dir=`dirname $f`

channelname=channel$((c/2 + 1))
mesondir=`echo ${name%%-dt*}`
dtdir=`echo $name | cut -d "-" -f 4 | cut -d "_" -f 1`

totalpath=./bag/$lat/$mesondir/$dtdir/
echo $name

mkdir -p $totalpath

cp $f $totalpath/$mesondir-$channelname.dat


done
done
done
