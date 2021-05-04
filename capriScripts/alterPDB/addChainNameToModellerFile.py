#!/bin/python
# -*- coding: utf-8 -*-
import os,sys,string

# Make sure the chain names match the ones in the PDB!

pdbFile=sys.argv[1]
chainName=sys.argv[2]
outputFile=sys.argv[3]


mdlFile=open(pdbFile,'r')
mdlLines=[line  for line in mdlFile.readlines()]
mdlFile.close()

of=open(outputFile,'w')
for ln in mdlLines:
     if ln.startswith("ATOM") or ln.startswith("TER"): # ATOM line
	   ln=ln[0:21]+chainName+ln[22:len(ln)] # Note that chain ID is the 22nd char in PDB sequences. Note also that :end means till char end-1
	   of.write(ln)  
 
     else: 
	   of.write(ln)


of.close()
