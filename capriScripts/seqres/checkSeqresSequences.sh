#!/bin/bash

seqresFile=$1
structureFile=$2

cut -c 20-70 $seqresFile > seqres_residues

python check_seqres_structure_seqs.py seqres_residues  $structureFile

