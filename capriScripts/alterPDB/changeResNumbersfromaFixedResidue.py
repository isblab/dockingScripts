#!/bin/python
# -*- coding: utf-8 -*-
import os,sys,string
pdbFile=sys.argv[1]
outputFile=sys.argv[2]
resTochangeFrom=int(sys.argv[3])

pdb=open(pdbFile,'r')
pdb_lines=[line.strip('\n') for line in pdb]
pdb.close()


of=open(outputFile,'w')
for ln in pdb_lines:
	ln_fields=ln.split()
	
	if len(ln_fields)<6 :
		of.write(ln+"\n")
		continue
	res=int(ln_fields[5])
	
	if res<resTochangeFrom:
	       of.write(ln+"\n")
	       continue
	else:
	      ln=ln.replace(ln_fields[5],str(res+1))
	      of.write(ln+"\n")
of.close()



