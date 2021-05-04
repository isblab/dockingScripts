#!/bin/bash

headerFile=$1
modelsList=$2
modelsHashTagFile=$3

for i in `cat $modelsList`
do
	mdlNum=`echo $i | sed s/m//`
	echo $mdlNum

	grep " "$mdlNum" MD5" $modelsHashTagFile >> $headerFile
done




