#include "system.h"
#include "ran250.h"
#include "preprocessor.h"
#include <math.h>
#include <time.h>
#include <float.h>

void Initialize(){
    // setup rounding mode

    int count=0;
    printf("********DENDRIMER SIM********\n");
//  Using
    printf("Using:\n");
#if USE_PBC==1
    printf("%1d. Periodic Boundary Conditions\n",count);
#endif
#if USE_CELL_LIST==1
    printf("%1d. Cell Lists\n",++count);
#endif
#if LJ==1
    printf("%1d. Morse - LJ approximation\n",++count);
#endif
    if ((count-1)<=0) printf("NO OPTIMIZATION\n");
    
    SetParams();
    Allocate();
    DendSetup();
    Dendrimer *pd=NULL;
    VecR vecrt={0.0, 0.0, 0.0};
//  INITIAL CONFIGURATION
//  RESTART IS NOT USED HERE
    for (int i=0;i<numOfDendrimers;i++){
   	    pd = dendrimer+i;
   	    DendDefineBondsTypes(pd);
   	}
    if (numOfIniConfigFiles==0) {
    	for (int i=0;i<numOfDendrimers;i++){
    	    pd = dendrimer+i;
			DendInitConfig(i,&vecrt); // calls Dendsetnpos 
       		RandomVector(&vecrt,DendMaxRad(pd));
    	}
	}else{
		DendReadCoords(fileIniConf,0); //calls Dendsetnpos or Dendsetpos
    }
    // Write Bonds File;
    FILE *pfBonds=fopen("bonds.dat","w");
    DendWriteBonds(pfBonds);
    fclose(pfBonds);
#ifdef LUT_MORSE
	LUTMorseOutputTables("morse.dat",0.0,4.0,0.01);
#endif
#ifdef LUT_FENE
	LUTFeneOutputTables("fene.dat",0.0,6.0,0.01);
#endif
}

void SetParams(void){
    if (seed==0){
        r250_init(time(NULL));
    }else if (seed!=0){
        r250_init(seed);
    }
    // G,F are defined in ReadIputData
    numOfMonomersInDend=2*(1-(int)pow((F-1),G+1))/(2-F);
    // eg.F=3
    // G2:14, G3:30, G4:62, G5:126, G6:254, G8:1022
    //
/*
    for (int ngen=0;ngen<=G;ngen++){
            numOfMonomersInDend +=2*(int) pow((F-1),ngen);
    }
*/
    numOfTotalMonomers  = numOfDendrimers*numOfMonomersInDend;
    NumofDispPerCycle   = numOfTotalMonomers;
    // 2000 frames for g6f3 dendrimers 16MB
    // 100.000 frames 16x50 = 800MB
    Density =(double)numOfTotalMonomers/(Box.x*Box.y*Box.z);
    RsqCut = (double*)malloc(N_TYPE_MORSE*sizeof(double));
    // maximum cutoff radius of all morse potentials
//#if USE_CELL_LIST==1
//    maxRc = MAX3(MrsCC.RCut,MrsCS.RCut,MrsSS.RCut);
//    maxRc = CellListSizeMaxrad * maxRc;
//    CellListSetup(Box, maxRc);
//#else
//    maxRc = -2.0;
//#endif
//    
    RsqCut[CC] = SQR(MrsCC.RCut); RsqCut[CS] = SQR(MrsCS.RCut);
    RsqCut[SS] = SQR(MrsSS.RCut);

    printf("numofmonomer in dendrimer:%-3d\n",numOfMonomersInDend);
    printf("Total number of monomers:%-4d\n",numOfTotalMonomers);
    printf("Density=%-lf\n\n",Density);

    Beta = 1.0 / Temp;
	printf("Temp=%lf,Beta=%lf\n",Temp,Beta);
}

void Allocate(){
/*
    CAUTION: matrices are allocated using G & F (generation
    Functionality) defined on readinput.c
    allocate dendrimer* Matrix
*/
    dendrimer=(Dendrimer*)malloc(numOfDendrimers*sizeof(Dendrimer));
    for (int dend=0;dend<numOfDendrimers;dend++){
        dendrimer[dend].numOfMonGen = (int*)malloc((G+1)*sizeof(int));
        dendrimer[dend].monomer     = (Monomer*)malloc(numOfMonomersInDend*sizeof(Monomer));
    }

/*
    allocate monomer* Matrix
    monomer: contains monomers of all dendrimers
*/
    for (int ndend=0;ndend<numOfDendrimers;ndend++){
        for (int mon=0;mon<numOfMonomersInDend;mon++){
            dendrimer[ndend].monomer[mon].bond           = (Monomer**)malloc(F*sizeof(Monomer*));
            dendrimer[ndend].monomer[mon].bondInterType  = (int*)malloc(F*sizeof(int));
        }
    }
	// Allocate LookupTables
#ifdef LUT_MORSE
	LUTMorseAllocate();
	// shift potential at cut
	LUTMorseBuild(TRUE);
#endif
#ifdef LUT_FENE
	LUTFeneAllocate();
	LUTFeneBuild(1.5);
#endif
}
