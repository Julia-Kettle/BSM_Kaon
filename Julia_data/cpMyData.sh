#!/bin/bash

for lat in 48cubed 64cubed
do


for f in ../../Fits/${lat}/mass/mass*/mass*_boots.dat
do
#echo $f
name=`basename $f`
dir=`dirname $f`

#echo $name
#echo $dir

newname=`echo $name | sed 's@_boots.dat@@g'`

totalpath=./mass/$lat/$newname/

mkdir -p $totalpath

cp $f $totalpath$newname.dat

done

for f in ../../Fits/${lat}/mass/mass*/decay*_boots.dat
do
#echo $f
name=`basename $f`
dir=`dirname $f`

#echo $name
#echo $dir

newname=`echo $name | sed 's@_boots.dat@@g'`

totalpath=./mass/$lat/$newname/

mkdir -p $totalpath

cp $f $totalpath$newname.dat

done

for c in 0 2 4 6 8
do
for f in ../../Fits/${lat}/bag/channel${c}/bag*/*_boots.dat
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
