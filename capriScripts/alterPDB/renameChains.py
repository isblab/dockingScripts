#!/bin/python
# -*- coding: utf-8 -*-

import os, sys,string

def renameChains():
     inputFile=sys.argv[1]

     inputChain=sys.argv[2]

     outputFile=sys.argv[3]

     outputChain=sys.argv[4]
	
     inp=open(inputFile,'r')
     inpLines=[line  for line in inp.readlines()]
     inp.close()

     out=open(outputFile,'w')

     for ln in inpLines:
          if ln.startswith("ATOM") or ln.startswith("TER"): # ATOM line

               ln=ln[0:21]+outputChain+ln[22:len(ln)] # Note that chain ID is the 22nd char in PDB sequences. Note also that :end means till char end-1        
               out.write(ln)

          else:
               out.write(ln)

     out.close()

# MAIN
if __name__ == '__main__':
    renameChains()


