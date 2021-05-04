#!/usr/bin/python

import os,sys,string

modelsFile=sys.argv[1]

outDir=sys.argv[2]

mcount=0

mf=open(modelsFile,'r')

for ln in mf:

	if ln.startswith("MODEL"):
		of=open(outDir+'/m'+str(mcount)+'.pdb','w')
		
		ln=mf.next()

		while not ln.startswith("ENDMDL"):
			of.write(ln)
			ln=mf.next()
	
		of.close()

		mcount=mcount+1
		print mcount

mf.close()














