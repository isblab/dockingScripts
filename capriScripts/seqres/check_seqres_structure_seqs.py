# -*- coding: utf-8 -*-
# This program uses the Biopython distribution. (http://www.biopython.org)

# 01/28/2010
# @Author Shruthi Viswanath

import os,sys,math,numpy,string
from subprocess import *
from optparse import OptionParser
import glob
from Bio.PDB import *

def getAA_Alphabet(resnam):
	
	if resnam=='ALA':
		return 'A'
	elif resnam=='CYS':
		return 'C'
	elif resnam=='ASP':
		return 'D'
	elif resnam=='GLU':
		return 'E'
	elif resnam=='PHE':
		return 'F'
	elif resnam=='GLY':
		return 'G'
	elif resnam=='HIS':
		return 'H'
	elif resnam=='ILE':
		return 'I'
	elif resnam=='LYS':
		return 'K'
	elif resnam=='LEU':
		return 'L'
	elif resnam=='MET':
		return 'M'
	elif resnam=='ASN':
		return 'N'
	elif resnam=='PRO':
		return 'P'
	elif resnam=='GLN':
		return 'Q'
	elif resnam=='ARG':
		return 'R'
	elif resnam=='SER':
		return 'S'
	elif resnam=='THR':
		return 'T'
	elif resnam=='VAL':
		return 'V'
	elif resnam=='TRP':
		return 'W'
	elif resnam=='TYR':
		return 'Y'
	else:
		return 'G'
	# default is GLYcine	
	
def makeSeqFile():
	
	seqresFile=sys.argv[1] # get the sequence part of seqres only, not the SEQRES AA_NUM etc, use Unix cut command. 
	
	structureFile=sys.argv[2] 
	
	seqres_seq=""
	seqres=open(seqresFile,'r')
	for line in seqres.readlines():
	     
	       residues_in_line=line.strip('\n').split(' ')
	       
	       for res in residues_in_line:
		    
		    if res=='':
			 continue
		    
		    alphabet=getAA_Alphabet(res)
		    seqres_seq=seqres_seq+alphabet
	seqres.close()
	
	# Get the structure
	parser=PDBParser()
	structure=parser.get_structure('pdbStruct', structureFile)
	
	# Get the list of residues in current structure
	structure_seq=""
	for res in structure.get_residues():
		nam=res.get_resname()
		alphabet=getAA_Alphabet(nam)
	        structure_seq=structure_seq+alphabet


	if seqres_seq!=structure_seq:
	       print "The SEQRES sequence and the sequence in the structure are not the same!!!"
	       print seqres_seq
	       print structure_seq

	else:
	       print "Both the sequences in SEQRES and the structure are the same. Phew!!"
	       print seqres_seq
	
    
# MAIN
if __name__ == '__main__':
    makeSeqFile()
