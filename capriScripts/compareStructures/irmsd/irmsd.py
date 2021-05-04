#!/usr/bin/python
# -*- coding: utf-8 -*-

# This program uses the Biopython distribution. (http://www.biopython.org)

# 11/30/2011
# @Author Shruthi Viswanath

import os,sys,math,numpy,string
import itertools
from Bio.PDB import *
from Bio import SeqIO,pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import warnings
from operator import itemgetter
       
''' 
I-rmsd as defined here: 
Mendez et al.
PROTEINS: Structure, Function and Genetics 52:51-67 (2003)
Assessment of Blind Predictions of Protein-Protein Interactions: Current Status of Docking Methods.
'''

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
		return ''
	# default is GLYcine	
	
def parseChainResiduesFromStructure(modelFile,chain):
     
        parser=PDBParser() # needed for parsing the structure
	ppb=PPBuilder() # needed for getting the sequence from a structure  
	
	modelStructure=parser.get_structure('pdbStruct',modelFile) # parse model structure
		
	resList=[]
	resSeq=""
		
	# Get the list of residues of the structure
	for res in modelStructure[0][chain]:
	       if res.has_id('CA') and res.has_id('N') and res.has_id('C') and res.has_id('O'):
		    resList.append(res) # Biopython syntax for getting the receptor chain. Structure->Model->Chain->Residue->Atoms
		    resSeq=resSeq+getAA_Alphabet(res.get_resname())
			    		    	       
	return(resList,resSeq)	 
	
def getSymmetricChainsList(modelFile,modelChains):
	
	symmetricChainsList=[] # e.g. [['A','B'],['B','A']] if A and B are identical chains
	symmetricChainsStrings=[] # e.g. [AB,BA]
	
	chainSequences={}
	chainsDict={}
		
	# Step 1. Get the sequences for all chains given and initialize dictionaries of symmetric chains
	for ch in modelChains:
	       (chResList,chSeq)=parseChainResiduesFromStructure(modelFile,ch)
	       
	       chainSequences[ch]=chSeq 
	
	       chainsDict[ch]=[ch]
	       
	# Step 2: Get identical chains in a dictionary. e.g. if chains A and B are identical in the string ABC, it returns {'A':['A','B'],'C':['C']}
	for indx1 in range(len(modelChains)):
	     ch1=modelChains[indx1]
	     
	     for indx2 in range(indx1+1,len(modelChains)):
		  ch2=modelChains[indx2]
		  
		  if ch1 in chainSequences and ch2 in chainSequences: # keys exist
		  
		       if almostIdentical(chainSequences[ch1],chainSequences[ch2]): # sequences identical
			      chainsDict[ch1].extend(chainsDict[ch2])
			      del chainsDict[ch2]			      
		       
	# Step 3: Initialize the strings of chains list
	numSymmChainStrings=0
	for ch in chainsDict:
	     if len(chainsDict[ch])>1:
		  numSymmChainStrings=numSymmChainStrings+math.factorial(len(chainsDict[ch]))
	
	if numSymmChainStrings==0: # no symmetric chains, return the chain list as it is 
	     return([modelChains])
	
	else:
	     for count in range(numSymmChainStrings):
		    symmetricChainsList.append(list(itertools.repeat(0,len(modelChains))))
	
	# Step 4: Get all the permutations of the new chain lists 
	for ch in chainsDict:
	     currSymmList=chainsDict[ch]
	
	     currIndicesList=getChainIndices(currSymmList,modelChains)
	     
	     sclcount=0
	     while sclcount<numSymmChainStrings: # loop over each string in the chains list. 
	     
		    for perm in itertools.permutations(currSymmList): # loop over each permutation of the current group, eg.(A,B) is one permutation, (B,A) is another
			 for i in range(len(currIndicesList)): # loop over the length of the string 
			      indx=currIndicesList[i]
			      symmetricChainsList[sclcount][indx]=perm[i]
	
			 sclcount=sclcount+1
	
	# Step 5: Convert the list to strings
	for chlist in symmetricChainsList:
	       chlist = ''.join(chlist)
	       symmetricChainsStrings.append(chlist)		
	
	return(symmetricChainsStrings)

def getChainIndices(currChainsList,actualChainString):
     
        return([actualChainString.index(ch) for ch in currChainsList])
        
def almostIdentical(seq1,seq2):
	matrix=matlist.blosum62  # Matrix for pairwise alignment
	pairAlignment=pairwise2.align.globalds(seq1,seq2, matrix, -10, -1) #gap open:-10, gap extend: -1
	
	numMismatches=len([indx for indx in range(len(pairAlignment[0][0])) if pairAlignment[0][0][indx]!=pairAlignment[0][1][indx]]) 
	
	if numMismatches<=6:
	     return True
	
	else:
	     return False
	     
def getCommonResiduesList(refFile,refChains,modelFile,modelChains):
        
        refComponent=[]
        modelComponent=[]
        
        matrix=matlist.blosum62  # Matrix for pairwise alignment
      
	for refch,mdlch in zip(refChains,modelChains): 
	       
	     (modelResidueList,modelSequence)=parseChainResiduesFromStructure(modelFile,mdlch)
	     (refResidueList,refSequence)=parseChainResiduesFromStructure(refFile,refch)	 
	  
	     modelRefAlignment=pairwise2.align.globalds(modelSequence, refSequence, matrix, -10, -1) #gap open:-10, gap extend: -1
	  	     	    	     
	     modelAligned=modelRefAlignment[0][0]
	     refAligned=modelRefAlignment[0][1]
	     	    
	     # this list holds all the residues that need to be deleted as they are missing in the other structure
	     modelResiduesToRemove=[] 
	     refResiduesToRemove=[]
	     
	     # this keeps track of the residue counts
	     mcount=0
	     rcount=0
	     
	     # iterate over the alignment looking for gaps. Both modelAligned and refAligned are of same length
	     for i in range(len(modelAligned)):
		  		       
		  if modelAligned[i]!='-' and refAligned[i]!='-':  # and modelAligned[i]==refAligned[i]:
		       rcount=rcount+1
		       mcount=mcount+1		       
		     
		  elif modelAligned[i]=='-':
		       refResiduesToRemove.append(refResidueList[rcount]) 
		       rcount=rcount+1
		       
		  elif refAligned[i]=='-':
		       modelResiduesToRemove.append(modelResidueList[mcount])
		       mcount=mcount+1
		         
	     # execute the removal in model and reference residue lists 
	     for itm in modelResiduesToRemove:
		  modelResidueList.remove(itm)
		 
	     for itm in refResiduesToRemove:
		  refResidueList.remove(itm)	
			  
	     # add in the lists to the dictionary represented by the component
	     refComponent.extend(refResidueList)
	     
	     modelComponent.extend(modelResidueList)	     
	  	     		   			    		    	       
	return(refComponent,modelComponent)
     
def getInterfaceResidues(refReceptorResList,refLigandResList,modelReceptorResList,modelLigandResList):
        ''' Returns a list of the interface residues of a given structure. '''
        
	interDis=10.0 # required interface distance between 2 atoms 
        
        safeDis=20.0  # the distance of CAs of 2 residues above which you can safely assume that no atoms of the 2 residues are within interDis #The actual distance is 21.4 A = 10+5.7*2(5.7 is the approximate distance from CA to farthest side-chain atom of TRP)
        
        refRecInterface={}
        refLigInterface={}
        
        modelRecInterface={}
        modelLigInterface={}
		
	# Find the interface residues	
				    
	for rrescount in range(len(refReceptorResList)):
	       rres=refReceptorResList[rrescount]
	       
	       for lrescount in range(len(refLigandResList)):
	       
		    lres=refLigandResList[lrescount]
					
		    if lres['CA']-rres['CA']>safeDis: # if the two residue CA's are over safeDis, then no atom in these 2 residues can be within interDis. 
		    	 continue
		    
		    elif lres['CA']-rres['CA']<interDis: # interface residues, since CA's are within 10 A
			 
			 if rrescount not in refRecInterface:
			      refRecInterface[rrescount]=rres['CA'].get_coord()	
			      modelRecInterface[rrescount]=modelReceptorResList[rrescount]['CA'].get_coord()
			   			   			     		      
			 if lrescount not in refLigInterface:
			      refLigInterface[lrescount]=lres['CA'].get_coord()	
			      modelLigInterface[lrescount]=modelLigandResList[lrescount]['CA'].get_coord()
						      
			 continue
		    
		    else: # calculate all atom distances only for residues' CA between 10 and 20 A. 
			 gotonext=False  # flag variable to get out of double for loop
			 for ratm in rres:
			      if not gotonext:
				   for latm in lres:
					if ratm-latm<interDis:					 
						  if rrescount not in refRecInterface:
						       refRecInterface[rrescount]=rres['CA'].get_coord()	
						       modelRecInterface[rrescount]=modelReceptorResList[rrescount]['CA'].get_coord()
																		       
						  if lrescount not in refLigInterface:
						       refLigInterface[lrescount]=lres['CA'].get_coord()	
						       modelLigInterface[lrescount]=modelLigandResList[lrescount]['CA'].get_coord()
						   						       
						  gotonext=True
						  break
			      else:
				   break	
	
		
	referenceInterface=combineCoords(refRecInterface,refLigInterface)	
	   
	modelInterface=combineCoords(modelRecInterface,modelLigInterface)

	return(referenceInterface,modelInterface)	
	
def combineCoords(receptorInterface,ligandInterface):
     
	coords=[]
	
	coords=[receptorInterface[indx] for indx in sorted(receptorInterface.iterkeys())]
		     
	coords.extend([ligandInterface[indx] for indx in sorted(ligandInterface.iterkeys())])
	
	coords=numpy.array(coords)
		    	     
	return(coords)         

def getRMSD(ra,rb):
        ''' ra and rb are NX3 arrays. '''
        
        R=numpy.zeros((3,3))
        b=numpy.zeros((3,3))
        
        lenSuperpose=ra.shape[0] # number of atoms 
        
        # Step 1. Translation
        centroidA=ra.sum(axis=0)/lenSuperpose # sum along the columns
        centroidB=rb.sum(axis=0)/lenSuperpose # sum along the columns 
       
        ra=ra-centroidA
        rb=rb-centroidB
              
        # Step 2. Get the matrix R
        for i in range(3):
	       for j in range(3):
		    R[i,j]=(rb[:,i]*ra[:,j]).sum()
		
        # Step 3. Get the eigen values of R, mu, and the eigen vectors a
        Rtranspose=R.conj().transpose()
        musquare=numpy.linalg.eig(numpy.dot(Rtranspose,R))[0]  # dot is the command for matrix multiply, this returns the eigen values in musquare
        a=numpy.linalg.eig(numpy.dot(Rtranspose,R))[1]
        
	mu=numpy.sqrt(musquare) 

	# Step 4: Get bk's
	for k in range(3):
	       b[:,k]=(numpy.dot(R,a[:,k]))/mu[k] 
	
	# Step 5: Get U
	[U,mu]=calcRotMatAndEigenValues(a,b,mu)

	# Step 6: Get distance D, and rmsd
        Dsquare=numpy.power(ra,2).sum()+numpy.power(rb,2).sum()-2.0*(mu.sum()) 

	if Dsquare > 0.0 :
		rmsd=math.sqrt(Dsquare)/math.sqrt(lenSuperpose) 
	else:
		rmsd=0.00           

	return(rmsd)

def calcRotMatAndEigenValues(a,b,mu):
	U=numpy.zeros((3,3))
		
	for i in range(3):
	      U=U+numpy.outer(b[:,i],a[:,i]) 
	     
	# Get the determinant of U
	detU=numpy.linalg.det(U)

	# check the determinant of U for reflection. 
	if numpy.abs(detU+1.0000)<0.001:  # Inversion detected 
	       minval,indx = mu.min(0),mu.argmin(0)
	       mu[indx]=-mu[indx]
	       
	       U=U-2*numpy.outer(b[:,indx],a[:,indx]) 
		
	return(U,mu)

    
def getIRMSD():
   
	warnings.filterwarnings('ignore')

   	modelFile=sys.argv[1] # model PDB file

	modelRch=sys.argv[2] # model receptor chain. Multiple chains are entered as a single string e.g AB

        modelLch=sys.argv[3] # model ligand chain. 
        
        refFile=sys.argv[4] # reference PDB file 
        
        refRch=sys.argv[5] # reference receptor chain
        
        refLch=sys.argv[6] # reference ligand chain
     
        # Step 0. Get the combinations of symmetric chains. CAUTION need to have calculated this beforehand and stored in a file. 
        if len(modelRch)>1:
	     mdlReceptorChainCombinations=getSymmetricChainsList(modelFile,modelRch)
	else:
	     mdlReceptorChainCombinations=[modelRch]
     
        if len(modelLch)>1:
	     mdlLigandChainCombinations=getSymmetricChainsList(modelFile,modelLch)        
        else:
	     mdlLigandChainCombinations=[modelLch]
	     
	least_irmsd=100000.0
	
	for mdlRchain in mdlReceptorChainCombinations:
	     
	    for mdlLchain in mdlLigandChainCombinations:
		 
	       # Step 1. Get the alignment between corresponding chains to eliminate residues missing in either file. 
	       (refReceptor,modelReceptor)=getCommonResiduesList(refFile,refRch,modelFile,mdlRchain)

	       (refLigand,modelLigand)=getCommonResiduesList(refFile,refLch,modelFile,mdlLchain)  

	       # Step 2. Get the interface residues according to the reference
	       (refCoords,modelCoords)=getInterfaceResidues(refReceptor,refLigand,modelReceptor,modelLigand) 
	     
	       if len(refCoords)==0 or len(modelCoords)==0:
			curr_irmsd=100.0
		
	       else:
			# Step 3. Compute RMSD		
	       		curr_irmsd=getRMSD(modelCoords,refCoords)
	       
	       if curr_irmsd<least_irmsd:
		    least_irmsd=curr_irmsd
        
        print "%.3f" %(least_irmsd)
              
# MAIN
if __name__ == '__main__':
    getIRMSD() 
