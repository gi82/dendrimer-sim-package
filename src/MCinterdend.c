/*
 * File:   MCinterdend.c
 * Author: georgiou
 *
 * Created on July 20, 2011, 5:25 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include "ran250.h" 
#include "system.h"
#include <math.h>
#include <ctype.h>
#include <float.h>

int main(int argc, char** argv) {
    ReadSimParams(argv[1]);
	ReadPotential(PotType,filePotential);
	OutputInputParams();

    Initialize();
//  if restart then DON'T TOUCH THE BOX SIZE
//  When both initail config & restart are used
//  system changes box size
    if (Restart==0){
		VecR vecr;
		Dendrimer *pd;
		// create large enough box to fit all dendrimers
		VecI tempBox ={110,110,110};
		VecI Boxbackup = {Box.x, Box.y, Box.z};
		// setting global Box to tempBox;
		Box.x = tempBox.x;	
		Box.y = tempBox.y;
		Box.z = tempBox.z;
		// generating dendrimers inside temp box
		printf("Generating dendrimers inside equilibration box\n");
		printf("Box size:%4d %4d %4d\n",Box.x,Box.y,Box.z);
		for (int dendid=0;dendid<numOfDendrimers;dendid++){
			pd     = dendrimer+dendid;
			vecr.x = (double)Box.x*dr250() - Box.x/2.0 ;
			vecr.y = (double)Box.y*dr250() - Box.y/2.0;
			vecr.z = (double)Box.z*dr250() - Box.z/2.0;
			// calculate position on the box
			WriteVecR(vecr,stdout);
			DendCalcCMass(pd);// uses npos
			DendPlaceCMass(pd,vecr);//uses npos
			WriteVecR(pd->CMass,stdout);
			DendSetpos(pd); //uses pos
		}
		fflush(stdout);
	
#if USE_PBC==1
		ApplyBoundaryConditions();
#endif
	
#if USE_CELL_LIST==1
	    maxRc = MAX3(MrsCC.RCut,MrsCS.RCut,MrsSS.RCut);
	    maxRc = CellListSizeMaxrad * maxRc;
	    CellListSetup(Box, maxRc);
		CellAlloc();
		CellNewBuild();
#else
	    maxRc = -2.0;
#endif
	
		// RUN EQUILIBRATION MONTE CARLO 
		MCequil(1,1000);		
		// RESTORE ORIGINAL BOX
		//
		Box.x = Boxbackup.x;
		Box.y = Boxbackup.y;
		Box.z = Boxbackup.z;
	
		printf("Moving CMass of dendrimers to normal Box\n");
		printf("Box size:%4d %4d %4d\n",Box.x,Box.y,Box.z);
	
		for (int dendid=0;dendid<numOfDendrimers;dendid++){
			pd     = dendrimer+dendid;
			vecr.x = (double)Box.x*dr250() - Box.x/2.0 ;
			vecr.y = (double)Box.y*dr250() - Box.y/2.0;
			vecr.z = (double)Box.z*dr250() - Box.z/2.0;
			printf("random pos->");
			WriteVecR(vecr,stdout);
			// no need to set npos to pos because of the
			// previous run of mc
			DendCalcCMass(pd);
			DendPlaceCMass(pd,vecr);
			printf("CMass  pos->");
			WriteVecR(pd->CMass,stdout);
			DendSetpos(pd);
		}
		fflush(stdout);
	}
	WriteBondVmd("bonds_vmd.dat");
#if USE_PBC==1
	ApplyBoundaryConditions();
#endif
#if USE_CELL_LIST==1
	CellDeAlloc();
	maxRc = MAX3(MrsCC.RCut,MrsCS.RCut,MrsSS.RCut);
	maxRc = CellListSizeMaxrad * maxRc;
	printf("Cell List range:%lf\n",maxRc);
	CellListSetup(Box, maxRc);
	CellAlloc();
	CellNewBuild();
#else
	maxRc = -2.0;
#endif
		// run final MC
	printf("energy:%lf\n",EnergySystem2());
	MC();   
	char topfname[100];	
	char coeffname[100];
	sprintf(topfname,"top.Dpbc_N%3d_%3d_D%3d.txt",numOfDendrimers,G,PotType);
	sprintf(coeffname,"coeff_N%d_%d_D%d.txt",numOfDendrimers,G,PotType);
	DendWriteinlmp(topfname,coeffname,1);
	sprintf(topfname,"top.Dnopbc_N%d_%d_D%d.txt",numOfDendrimers,G,PotType);
	DendWriteinlmp(topfname,coeffname,0);

    return (EXIT_SUCCESS);
}
