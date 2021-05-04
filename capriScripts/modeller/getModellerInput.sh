#!/bin/bash

tgt=$1
template=$2
aliFile=$3 # the part from the brk file which contains the alignment, REMARKS can be included at the beginning of the line
outAli=$4  

template_pdb=`echo $template | cut -d '_' -f 1`
template_chain=`echo $template | cut -d '_' -f 2`

# STEP 1 get sequence file 
echo ">"$tgt > $tgt.seq
grep EMBOSS $aliFile| awk '{print $4}' >> $tgt.seq

# STEP 2 get alignment file
tgtSeq=`grep EMBOSS $aliFile|cut -d ' ' -f 4`
templateSeq=`grep SWALL $aliFile | cut -d ' ' -f 4` 

first=`grep ATOM $template.pdb | awk '{print $6}' | head -n 1`
last=`grep ATOM $template.pdb | awk '{print $6}' | tail -n 1`

echo "C; "$tgt" -> "$template > $outAli
echo ">P1;"$template_pdb >> $outAli
printf "structureX:"$template".pdb:$first: :$last: : : : :" >> $outAli

for ln in `echo $tgtSeq`
do 
	printf "\n%s" $ln >> $outAli
done

echo "*" >> $outAli

echo ">P1;"$tgt >> $outAli
printf "sequence:"$tgt": : : : : : : :" >> $outAli

for ln in `echo $templateSeq`
do 
	printf "\n%s" $ln >> $outAli
done

echo  "*" >> $outAli

# STEP 3 get Modeller script
echo "from modeller import *" > modeller_input.py 
echo "from modeller.automodel import *" >> modeller_input.py 
echo "env = environ()" >> modeller_input.py 
echo "a = automodel(env,alnfile  = '"$outAli"',knowns='"$template_pdb"',sequence = '"$tgt"')" >> modeller_input.py 
echo "a.starting_model= 1" >> modeller_input.py 
echo "a.ending_model  = 1" >> modeller_input.py 
echo "a.make()" >> modeller_input.py 

# STEP 4. Run Modeller

~/bin/modeller9v8/bin/mod9v8 modeller_input.py

# STEP 5. Rename modeller output

mv $tgt.B99990001.pdb $tgt.pdb

