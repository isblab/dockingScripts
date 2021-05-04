#!/bin/bash
pdb=$1
chain=$2

# get the atom number for the new OXT atom
lastAtomNum=`grep ^.....................$chain $pdb | tail -n 1 | awk '{print $2}'`  
nextAtomNum=`expr $lastAtomNum + 1`

# check for the last O atom, OXT is going to be close to that. So add a new OXT line with the co-ordinates equal to that of the last O atom
lastOxyAtom=`grep ^.....................$chain $pdb | grep 'ATOM.........O\ ' | tail -n 1` 

oxyAtomNum=`echo $lastOxyAtom | awk '{print $2}'`


# if incrementing the number of atoms changes the number of digits, e.g. 99 to 100
if [ ${#oxyAtomNum} != ${#nextAtomNum} ]
then 
	replaceAtomNumStr='.'$oxyAtomNum
else
	replaceAtomNumStr=$oxyAtomNum
fi

# make sure the PDB columns are properly aligned: if using echo all the spaces tend to get compressed to one space
grep ^.....................$chain $pdb | grep 'ATOM.........O\ '| tail -n 1 | sed -e s/$replaceAtomNumStr/$nextAtomNum/g -e s/\ O\ \ /\ OXT/g >> $pdb 

