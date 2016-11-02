/* 
 * File:   singledend.c
 * Author: giannis
 *
 * Created on August 2, 2011, 11:22 AM
 */

#include <stdio.h>
#include <stdlib.h>
#include "preprocessor.h"
#include "system.h"
#include <float.h>

int main(int argc, char** argv) {
    ReadSimParams(argv[1]);
	ReadPotential(PotType,filePotential);
	OutputInputParams();
    Initialize();
	if (Restart == 0){
    printf("max radius of dendrimer: %lf\n",DendMaxRad(dendrimer));
    VecR       vecr={0.0,0.0,0.0};
    Dendrimer *pd =dendrimer;
    
    DendSetnpos(pd);
    DendCalcCMass(pd);
    DendPlaceCMass(pd,vecr);
    WriteVecR(pd->CMass,stdout);
    DendSetpos(pd);
	}	
	FILE *pf;
	pf=fopen("totalpot.dat","w");
	fprintf(pf,"%13.6f %13.6f %13.6f\n",MrsCC.RCut,MrsCS.RCut,MrsSS.RCut);
	for (double rr=0.5;rr<=1.8;rr+=0.01){
		fprintf(pf,"%13.6f %13.6f %13.6f %13.6f\n",rr,Fene(&fnCC,rr),Morse(&MrsCC,rr),Fene(&fnCC,rr)+Morse(&MrsCC,rr));
	}
	fclose(pf);
#if USE_CELL_LIST==1
	    maxRc = MAX3(MrsCC.RCut,MrsCS.RCut,MrsSS.RCut);
	    maxRc = CellListSizeMaxrad * maxRc;
	    CellListSetup(Box, maxRc);
		CellAlloc();
		CellNewBuild();
#endif
    //printf("size=%ld, dig=%d, mant=%d, e=%lf", sizeof(double), DBL_DIG, DBL_MANT_DIG, DBL_EPSILON);
    MC();
    return (EXIT_SUCCESS);
}




