#!/bin/bash

modelsDir=$1

headerFile=$2

modelsPrefixList=$3

outFile=$4

preOutFile=$outFile'.pre'
rm $outFile 
# remove previous versions

cat $headerFile > $preOutFile

if [ $modelsPrefixList = 'NA' ]
then 
	numModels=10
	for mdlNum in `seq 0 9`
	do 	
		submittedModelNum=`expr $mdlNum + 1`	 
		echo "MODEL        "$submittedModelNum >> $preOutFile
		cat $modelsDir'/'m$mdlNum'.pdb' >> $preOutFile
		echo "ENDMDL"    >> $preOutFile
        done

else
	numModels=`wc -l $modelsPrefixList|awk '{print $1}'`
	mdlCount=0         

	for mdlNum in `cat $modelsPrefixList`
        do  
		submittedModelNum=`expr $mdlCount + 1`
		mdlCount=`expr $mdlCount + 1`	 
		echo "MODEL        "$submittedModelNum >> $preOutFile
		cat $modelsDir'/m'$mdlNum'.pdb' >> $preOutFile
		echo "ENDMDL" >> $preOutFile
	done
fi

# post processing
grep -v -e ^$ -e "REMARK   6" -e LOOPP -e "REMARK MODEL" -e "REMARK PARENT" -e "REMARK SCORE" -e "REMARK prediction" -e END$ $preOutFile > $outFile

echo "END" >> $outFile

echo "NOT done yet!!! Write TER for the tail added models manually!!!!"



