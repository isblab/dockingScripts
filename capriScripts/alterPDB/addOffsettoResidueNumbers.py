#!/usr/bin/python
# -*- coding: utf-8 -*-

# 01/28/2010
# @Author Shruthi Viswanath

import os,sys,math,string
	
def getPDB():
   
	  inputFile=sys.argv[1]
	  outputFile=sys.argv[2]
	  chain=sys.argv[3]
	  offset=int(sys.argv[4]) # by how many integers to move each residue
	  
	  inp=open(inputFile,'r')
	  pdbLines=[line for line in inp.readlines() ] # only ATOM lines are taken
	  inp.close()
	  
	  out=open(outputFile,'w')
	  
	  # write output file
	  for line in pdbLines:
		    if line.startswith('ATOM') and line[21]==chain:
			 old_res_number=line[23:27]
			 
			 new_res_number=str(int(old_res_number)+offset)
			 
			 if len(new_res_number)==1:
			      line=line[0:23]+"  "+new_res_number+line[26:len(line)]
			 elif len(new_res_number)==2:
			      line=line[0:23]+" "+new_res_number+line[26:len(line)]
			 else:
			      line=line[0:23]+new_res_number+line[26:len(line)]
			 			 
		         out.write(line)
		    

	  out.close()
				
# MAIN
if __name__ == '__main__':
    getPDB()
 
