#!/usr/bin/python
# -*- coding: utf-8 -*-

# This program uses the Biopython distribution. (http://www.biopython.org)

# 11/30/2011
# @Author Shruthi Viswanath

import os,sys,math,numpy,string
from Bio.PDB import *
from Bio import SeqIO
import warnings
from operator import itemgetter
import csv

''' 
Fraction of native contacts as defined here: 
Mendez et al 
PROTEINS: Structure, Function and Bioinformatics 60:150-169 (2005)
Assessment of CAPRI predictions in rounds 3-5: shows progress in docking procedures.
'''

def parseResiduesFromStructure(modelFile,chains):
     
        parser=PDBParser() # needed for parsing the structure
	ppb=PPBuilder() # needed for getting the sequence from a structure  
	
	modelStructure=parser.get_structure('pdbStruct',modelFile) # parse model structure
	
	chainList=[ch for ch in chains]

	resList=[]
		
	# Get the list of residues of the structure
	for chain in chainList:
	     for res in modelStructure[0][chain]:
		  if res.has_id('CA'):
		       resList.append(res) # Biopython syntax for getting the receptor chain. Structure->Model->Chain->Residue->Atoms
	
	return(resList)	
     
def getContactResidues(pdbFile,rChains,lChains):
        ''' Returns a list of the interface residues of a given structure. '''
        
	interDis=5.0 # required interface distance between 2 atoms 
        
        safeDis=20.0  # the distance of CAs of 2 residues above which you can safely assume that no atoms of the 2 residues are within interDis #The actual distance is 21.4 A = 10+5.7*2(5.7 is the approximate distance from CA to farthest side-chain atom of TRP)
		
	receptorResList=parseResiduesFromStructure(pdbFile,rChains)
	ligandResList=parseResiduesFromStructure(pdbFile,lChains)
	
	# Initialize the contact residue list
	contactResidues=[] # list of tuples of the type, (receptor_residue_number,receptor_residue_name,ligand_residue_number,ligand_residue_name)
	
	# Find the interface residues
	for rrescount in range(len(receptorResList)):
	    rres=receptorResList[rrescount]
	    
	    for lrescount in range(len(ligandResList)):
		  lres=ligandResList[lrescount]
		  
		  if lres['CA']-rres['CA']>safeDis: # if the two residue CA's are over safeDis, then no atom in these 2 residues can be within interDis. 
		      continue
		  
		  elif lres['CA']-rres['CA']<interDis: # interface residues, since CA's are within interDis
		       contactResidues.append((rres.get_resname(),rrescount,lres.get_resname(),lrescount)) # this way numbers are standardized! 
		       continue
		  
		  else: # calculate all atom distances only for residues' CA between 10 and 20 A. 
		       gotonext=False  # flag variable to get out of double for loop
		       for ratm in rres:
		 	    if not gotonext:
		 	      for latm in lres:
				   if ratm-latm<interDis:
					     contactResidues.append((rres.get_resname(),rrescount,lres.get_resname(),lrescount))
					     gotonext=True
					     break
			    else:
				 break
		  	
	return(contactResidues)	

def loadPIEMatrix(pieCSVFile):
         ''' Load PIE CSV file into a matrix. Takes the contact values from a file. This file was obtained from the table in DOCK/PIE JCP 2011 paper. '''
         
         mat=csv.reader(open(pieCSVFile,'r'),delimiter='\t')
         subsTable={}
     
         rowNum=0
         
         for row in mat:
	      if rowNum==0:
		   colIndices=row # starts with ' '
		   rowNum+=1
		   continue
	      
	      colNum=0
	      
	      for col in row:
		   if colNum==0:
			rowIndex=col
			colNum+=1
			continue
		   
		   subsTable[rowIndex,colIndices[colNum]]=col
		   colNum+=1
		   
	      rowNum+=1
	      
	 return subsTable 
         
def fnatContacts():
   
	warnings.filterwarnings('ignore')

   	modelFile=sys.argv[1] # model PDB file

	modelRch=sys.argv[2] # model receptor chain. Multiple chains are entered as a single string e.g AB

        modelLch=sys.argv[3] # model ligand chain. 
        
        pieCSVFile=sys.argv[4] # matrix which contains the amino acid- substitution values
        
        modelContactResidues=getContactResidues(modelFile,modelRch,modelLch)
        
        pieMatrix=loadPIEMatrix(pieCSVFile)
 
        #print len(modelContactResidues) 
	for (rres,rresc,lres,lresc) in modelContactResidues:
	        
		print rres,rresc,lres,lresc,pieMatrix[rres,lres]
		
		#TODO the output is not sorted. 
		#TODO Also the value of the contact energy is stored as a string not float. 

# MAIN
if __name__ == '__main__':
    fnatContacts()
 
