# -*- coding: utf-8 -*-

#!/usr/bin/python 

import os,sys,string,numpy

def fileWrite(dat,outFile):
     
	f=open(outFile,'w')
	for val in dat:
	  print >>f, str(val)
        f.close()
        
def fileRead(dataFile):
	dat=[]
	
	f=open(dataFile,'r')
        for ln in f.readlines():
	     dat.append(float(ln.strip()))
	  
        f.close()
        
        return dat

def rerankAtmPot():
        tgtName=sys.argv[1]
	potType=sys.argv[2] # e.g. atmver19.maxlhs 
	sortType=sys.argv[3] # pie and its variants use descending sort
		
	energyFile=tgtName+'.'+potType+'.energy' 
	outFile=tgtName+'.'+potType+'.indices.reranked.txt'
	
	# get energy values
	energy=numpy.array(fileRead(energyFile)) 
	
	# sort them, in ascending order: lower energy at first
	if sortType=="asc":
		ind=numpy.argsort(energy,axis=0)
	elif sortType=="desc":
		ind=numpy.argsort(energy,axis=0)[::-1]
	
	# model numbers not indices
	fileWrite(ind,outFile)
	
if __name__=='__main__':
	rerankAtmPot() 
