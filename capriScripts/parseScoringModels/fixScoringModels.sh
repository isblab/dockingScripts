#!/bin/bash

inFile=$1
outFile=$2

rm $outFile

# make residue numbers positive for all chains
for ch in A B
do 
       leastNegResNum=`grep ^ATOM.................$ch $inFile | head -1 | awk '{print $6}'`
       echo $leastNegResNum

       if [ $leastNegResNum -gt 0 ] 
       then 
		offset=0 
       else
      		 offset=`expr $leastNegResNum \* -1`
       		offset=`expr $offset + 1`
	fi

	python ~/capriScripts/alterPDB/addOffsettoResidueNumbers.py $inFile $inFile'.'$ch'.new' $ch $offset          
     
       cat $inFile'.'$ch'.new' >> $outFile	 		

	rm $inFile'.'$ch'.new'

done 

# change atom names!
sed -i -e s/OT1/O\ \ /g -e s/OT2/OXT/g $outFile

sed -i -e s/O1/O\ /g -e s/O2\ /OXT/g $outFile
