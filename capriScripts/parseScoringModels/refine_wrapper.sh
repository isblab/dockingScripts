#!/bin/bash

jobno='T59'
rchain='A'
lchain='B'

rm -r refine

mkdir refine

for i in `seq 0 1681`
do 
	echo "m"$i >> refine/model_list_all
done

cd refine

split -a 2 -d -l 20 model_list_all model_list_



for count in `cat ../seqnumbers_list`
do
     echo "#$ -N refine_$count
#$ -cwd
#$ -V
#$ -o refine_$count'.out' 
#$ -e refine_$count'.err' 
#$ -l h_rt=50:00:00
#$ -q all.q
#$ -pe mpich 1
#$ -S /bin/bash
set DISPLAY='' 
for mdl in \`cat model_list_$count\`
do 
	sh ~/data/capri28/scripts/fixScoringModels.sh ~/data/capri28/score/models/\$mdl'.pdb' ~/data/capri28/score/refine/\$mdl'.pdb'
done
sh ~/dockbypierr/processAndMinimizeModels.sh model_list_$count $jobno $rchain $lchain 

" > q_refine_$count
qsub q_refine_$count

done

