/*
Fast Artificial Neural Network Library (fann)
Copyright (C) 2003-2012 Steffen Nissen (sn@leenissen.dk)

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "fann-source/src/floatfann.c"
#include "fann-source/src/include/fann_cpp.h"

#include <ios>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <vector>
#include <fstream>
#include <algorithm>
using namespace std;

#define FILENAME_LEN 500
#define MAX_NUM_ATOMS 50000
#define MAX_LINE_LEN 200
#define FIELD_LEN 30

#define MAX_HYD_DONOR_COVALENT_BOND_LEN_SQUARE 1.2544 //1.12*1.12 

#define MAX_DISTANCE 8.0 
#define	MIN_DISTANCE 0.0 

#define NUM_PARTICLE_TYPES 4
	
#define NUM_BINS 2

#define NUM_PARAMS NUM_PARTICLE_TYPES*NUM_PARTICLE_TYPES*NUM_BINS
		
const float r1[NUM_BINS]={0.0,4.0}; // yes there is a hole in the middle
const float r2[NUM_BINS]={4.0,8.0};

#define FEATURE_LOWER_VAL -1.0
#define FEATURE_UPPER_VAL 1.0

#define NUM_INPUTS 32
#define NUM_OUTPUT 1


typedef struct atomInfo
{
	char atomName[FIELD_LEN];
	char residueName[FIELD_LEN]; // this should ideally be in residue structure but we are not storing these as classes so we cant access it! 
	float coord[3];
	struct atomInfo* donor;  // electro-negative atom bonded covalently to the atom. This field is NON-NULL only for hydrogens
} atom;	

typedef struct residueInfo // residue is just a set of atoms for us. 
{	vector<atom> atomList;
} residue;

// Note, this works only with MOIL PDB (polar hydrogen files) and HAAD files so far. 
bool isHydrogen(char *atom_name)
{
	if(atom_name[0] == 'H') 
		return true;
	else if(atom_name[0]>='1' && atom_name[0]<='9' && atom_name[1]=='H')
		return true;
	
	return false;
}

bool isElectroNegative(char *atom_name)
{
	if(atom_name[0] == 'O' || atom_name[0] == 'S' || atom_name[0] == 'N')
		return true;

	return false;
}

bool charBelongstoString(char str[],char c)
{
	for(int i=0;i<strlen(str);i++)
	{
		if(str[i]==c)
			return true;
	}
	return false;  

} 

float distanceSquare(float a[3], float b[3])
{
	return ((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));
}

class Molecule { // refers 
public:
	vector<atom> polarHyd; // list of polar hydrogens in the molecule
	vector<atom> electroNeg; // list of electronegative atoms in the molecule
	
	vector<residue> residueList;
	
	Molecule(string pdbFile,char chains[FIELD_LEN]);
	
	void createHBondAtomsList();
	
	void getHbondContacts(vector<atom> acceptorList,float cutoffDistance, float *numContacts);
	
	~Molecule() {};
};

Molecule::Molecule(string pdbFile,char chains[FIELD_LEN])
{
	// variables for various PDB fields
	char strbuf[MAX_LINE_LEN],line[MAX_LINE_LEN],ch;
	char atomCheck[FIELD_LEN],atmName[FIELD_LEN],resName[FIELD_LEN],xcoord[FIELD_LEN],ycoord[FIELD_LEN],zcoord[FIELD_LEN],resNumString[FIELD_LEN];

	float xpos,ypos,zpos;
	int i,j,resNum;
	
	ifstream fin(pdbFile.c_str()); // pdb file to be parsed
	
	int currResidueNum=9999; // unreachable residue

	// Parse the PDB file of the given model
	while(fin.getline(strbuf,MAX_LINE_LEN))
	{
			strcpy(line,strbuf); // all splicing operations are done on the string, line

			strncpy(atomCheck,line,4);
			atomCheck[4]='\0';

			if(strcmp(atomCheck,"ATOM")==0) // if a valid atom line, located an atom object
			{
				// get atom name
				for(i=0,j=12;j<16;j++)
				{	if(line[j]!=' ')
					{  atmName[i]=line[j];  
					   i++;
					}
				}
				atmName[i]='\0';
				
				if(!isHydrogen(atmName) && !isElectroNegative(atmName))
					continue;

				// get residue name
				for(i=0,j=17;j<20;i++,j++)
					resName[i]=line[j];
				resName[i]='\0';

				// get chain
				ch=line[21];
				
				if(!charBelongstoString(chains,ch)) // we are calling this function twice per molecule. Can call just once but this is neater. 
					continue;

				/* get residue number */
				for(i=0,j=22;j<26;i++,j++)
					resNumString[i]=line[j];
				resNumString[i]='\0';

				resNum=atoi(resNumString);
				
				if(currResidueNum!=resNum) // new residue begins here
				{
					currResidueNum=resNum;
					
					residue currResidue;
					residueList.push_back(currResidue); // push an empty residue to the residue list for the current molecule. 
				}	
							
				// get the X,Y,Z coordinates
				for(i=0,j=30;j<38;i++,j++)
					xcoord[i]=line[j];
				xcoord[i]='\0';
				xpos=atof(xcoord);

				for(i=0,j=38;j<46;i++,j++)
					ycoord[i]=line[j];
				ycoord[i]='\0';
				ypos=atof(ycoord);

				for(i=0,j=46;j<54;i++,j++)
					zcoord[i]=line[j];
				zcoord[i]='\0';
				zpos=atof(zcoord);

				// add atom to the current residue
				atom newatm;
				strcpy(newatm.atomName,atmName);
				strcpy(newatm.residueName,resName);
				newatm.coord[0]=xpos;
				newatm.coord[1]=ypos;
				newatm.coord[2]=zpos;
				newatm.donor=NULL; 
				
	
				// Add the current atom object to the vector of atom objects in the molecule object
				residueList.back().atomList.push_back(newatm);
										
			} //end ATOM line

	}  //end parsing in while loop
		
	fin.close();
}

void Molecule::createHBondAtomsList()
{
	
	for(vector<residue>::iterator rs=residueList.begin(); rs!=residueList.end();++rs) // for each residue in the molecule
	{
		for(int i=0;i<(*rs).atomList.size();i++) // for each atom in the residue
		{
			if(isHydrogen((*rs).atomList[i].atomName)) // find if the hydrogen is bonded to a polar atom. If so add the polar atom's info to the hydrogen atom and store the hydrogen in the polar hydrogen list
			{
				for(int j=0;j<(*rs).atomList.size();j++)
				{
					if(isElectroNegative((*rs).atomList[j].atomName)) // skip if examining the same atom or if the atom being examined is non-electronegative
					{
						if(distanceSquare((*rs).atomList[i].coord,(*rs).atomList[j].coord)<=MAX_HYD_DONOR_COVALENT_BOND_LEN_SQUARE ) // then it means the hydrogen is bonded covalently to the electronegative atom
						{
							(*rs).atomList[i].donor=&((*rs).atomList[j]);
							polarHyd.push_back((*rs).atomList[i]);
							
							// cout <<  (*rs).atomList[i].atomName << " " << (*rs).atomList[i].residueName << " " <<  (*rs).atomList[j].atomName << " " << (*rs).atomList[j].residueName << endl ;
							
							break; // we are done with this hydrogen atom!
						}
					
					}
					
				}
											
			}
				
		
			else if(isElectroNegative((*rs).atomList[i].atomName)) // else if the atom is an electronegative atom, add it to the list of electronegative atoms for the molecule
				electroNeg.push_back((*rs).atomList[i]);
			
		}
		
	} 
	
}

int getDistanceBin(float dist)
{       int i;
	
	if(dist>=MAX_DISTANCE||dist<MIN_DISTANCE)  //erroneous values
		return -1;
	
	for(i=NUM_BINS-1;i>=0;i--)
		if(dist>=r1[i] && dist<r2[i]) // going backwards because the lower bins are less populated
				return(i);
				
	return(-1);		   
		
}

int getParticleType(char resName[FIELD_LEN])
{
	// hydrophobic
	if(strcmp(resName,"ALA")==0 || strcmp(resName,"GLY")==0 || strcmp(resName,"VAL")==0 || strcmp(resName,"ILE")==0 || strcmp(resName,"LEU")==0 || strcmp(resName,"PRO")==0 || strcmp(resName,"TRP")==0 || strcmp(resName,"PHE")==0 || strcmp(resName,"MET")==0)
		return 0;

	//polar
	else if(strcmp(resName,"SER")==0 || strcmp(resName,"THR")==0 || strcmp(resName,"CYS")==0 || strcmp(resName,"ASN")==0 || strcmp(resName,"GLN")==0 || strcmp(resName,"TYR")==0)
		return 1;

	// positive
	else if(strcmp(resName,"ARG")==0 || strcmp(resName,"LYS")==0 || strcmp(resName,"HIS")==0 )
		return 2;
	
	// negative
	else if(strcmp(resName,"ASP")==0 || strcmp(resName,"GLU")==0 )
		return 3; 
	
	return -1;

}

void Molecule::getHbondContacts(vector<atom> acceptor, float cutoff,float *numContacts)
{
	for(vector<atom>::iterator hyd=polarHyd.begin();hyd!=polarHyd.end();hyd++)
	{
		for(vector<atom>::iterator acc=acceptor.begin();acc!=acceptor.end();acc++)
		{
			float hydbondLength=sqrt(distanceSquare((*hyd).coord,(*acc).coord));
			
			if(hydbondLength<cutoff)
			{
				
				int r=getDistanceBin(hydbondLength);

				int htype=getParticleType((*hyd).residueName);

				int atype=getParticleType((*acc).residueName);
		
				if(r >-1 && atype >-1 && htype >-1 )
				{
					int indx = NUM_BINS * (NUM_PARTICLE_TYPES * htype + atype) +r;   // map the particle and distance bin indices to the parameter index. 
			
					// cout << indx << " " << htype << " " << atype <<" " << r << endl ;
					numContacts[indx]+=1.0;
				}
			}
		}
	}
	
}

/* Scale the number of contacts */
void scaleFeatures(float *numContacts)
{
	// Hard-coded feature values
	float featureMin=0;
	float featureMax[NUM_INPUTS]={62,519,32,380,23,206,28,163,48,691,82,881,36,343,56,320,72,590,91,577,115,615,162,598,9,80,17,98,15,117,13,89};
	

	for(int i=0;i<NUM_INPUTS;i++)
	{			      
       		if(numContacts[i] == featureMin)
                	numContacts[i] = FEATURE_LOWER_VAL;
       		else if(numContacts[i] == featureMax[i])
                	numContacts[i] = FEATURE_UPPER_VAL;
        	else
                	numContacts[i] = FEATURE_LOWER_VAL + (FEATURE_UPPER_VAL-FEATURE_LOWER_VAL) * (numContacts[i]-featureMin)/(featureMax[i]-featureMin);
	}

}

/* Predict the hydrogen bond energy given the scaled interface hb features */
float nnEnergy(string net_file, float *scaled_hbfeatures)
{
	unsigned int i;

	fann_type *output;
	
	/* Step 1. Create network here */
	FANN::neural_net ann; 
	ann.create_from_file(net_file);
	
	/* Step 2. Get the test inputs loaded */
	fann_type* test_input = (fann_type *) calloc(NUM_INPUTS, sizeof(fann_type));

	for(i=0;i<32;i++)
		test_input[i]=scaled_hbfeatures[i];

	ann.reset_MSE(); // do this every time

        output=ann.run(test_input);

        return(output[0]);

}

int main(int argc,char *argv[])
{
	string pdbFile, netFile;
	char pdbReceptorChains[FIELD_LEN],pdbLigandChains[FIELD_LEN];

	pdbFile = argv[1];
	strcpy(pdbReceptorChains,argv[2]);
	strcpy(pdbLigandChains,argv[3]);
	netFile = argv[4];	

	float* numContacts=new float[NUM_INPUTS]; 
	
	// STEP 1a. Get all relevant atoms i.e. hydrogens and all electro-ve atoms in the appropriate residue lists 
	Molecule receptor(pdbFile,pdbReceptorChains);
	
	Molecule ligand(pdbFile,pdbLigandChains);

	// STEP 1b.  Create the list of polar hydrogen and electronegative atoms for the molecule  
	receptor.createHBondAtomsList();
	
	ligand.createHBondAtomsList();
	
	// STEP 2. Get the interface hydrogen-acceptor contact numbers by counting acceptors from both proteins. 
	for(int i=0;i<NUM_INPUTS;i++)
		numContacts[i]=0.0;
	
	receptor.getHbondContacts(ligand.electroNeg,MAX_DISTANCE,numContacts); 
	
	ligand.getHbondContacts(receptor.electroNeg,MAX_DISTANCE,numContacts); 

	// STEP 3. Now numContacts has the unscaled number of hydrogen bond contacts. Need to scale it
	scaleFeatures(numContacts);
	
	// STEP 4. Get the NN energy
	printf("%.4f\n", nnEnergy(netFile, numContacts));

	return 0;
}
