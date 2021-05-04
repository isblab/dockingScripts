#!/bin/bash

seqNumbersFile=$1
taskType=$2
target=$3

case $taskType in 
"piepisa")
		for i in `cat $seqNumbersFile`
		do 
			for scr in pie pisa
			do
				sed s/SeqNumber/$i/g ~/data/capri28/scripts/q_$scr > q_$scr'_'$i
				qsub q_$scr'_'$i
			done
		done ;;
"cx") 	sed s/tgt/$target/g ~/data/capri28/scripts/q_cx > q_cx
	qsub q_cx ;;

esac


