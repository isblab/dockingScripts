#!/usr/bin/python

import os,sys,string

inPdb=sys.argv[1]
outPdb=sys.argv[2]

residueList=['ALA','CYS','ASP','ARG','TRP','PHE','LEU','ILE','GLY','SER','THR','TYR','PRO','MET','HIS','GLN','ASN','GLN','VAL','LYS','GLU']

inf=open(inPdb,'r')
inputLines=[ln.strip() for ln in inf.readlines()]
inf.close()

outf=open(outPdb,'w')

for ln in inputLines:
	if ln.startswith('ATOM') and not (ln[17:20] in residueList):
		print ln[17:20]
		continue

	print >>outf, ln

outf.close()
