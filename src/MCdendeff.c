/* 
 * File:   main.c
 * Author: georgiou
 *
 * Created on April 18, 2011, 4:02 PMsrc_1
 */

#include <stdio.h>
#include <stdlib.h>
#include "ran250.h"
#include "system.h"
#include <math.h>
#include <ctype.h>
#include <float.h>
#include "preprocessor.h"

int main(int argc, char** argv) {
    ReadSimParams(argv[1]);
	ReadPotential(PotType,filePotential);
	OutputInputParams();
    Initialize();
    if (numOfDendrimers==2){
        printf("Simulating two Dendrimers\n");
        double maxCDdist, dendDist;
        VecR   dr;
        Dendrimer *pd1=NULL;
        Dendrimer *pd2=NULL;

        pd1 = dendrimer;
        pd2 = dendrimer+1;
        //output initial configurations
        DendWriteProps(pd1,stdout);
        DendWriteProps(pd2,stdout);
        // Calculate initial distance of CM of Dendrimers
        DendSetnpos(pd1);
        DendCalcCMass(pd1);
        DendSetnpos(pd2);
        DendCalcCMass(pd2);
        printf("Position dendrimer 1:"); WriteVecR(pd1->CMass,stdout);
        printf("Position dendrimer 2:"); WriteVecR(pd2->CMass,stdout);
        
        dr.x = pd1->CMass.x - pd2->CMass.x;
        dr.y = pd1->CMass.y - pd2->CMass.y;
        dr.z = pd1->CMass.z - pd2->CMass.z;
        
        dendDist  = sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.y));
        printf("initial distance:%8.3f\n",dendDist);

        // Calculate maximum allowed distance between monomers
        maxCDdist = DendMaxRad(pd2);
        printf("Max allowed distance:%8.3f\n",maxCDdist);
        if (iniDendDist>maxCDdist){
            printf("!!!!Initial input distance exceeds maximum distance :%8.3f\n",iniDendDist);
        }
        
        dr.x = -iniDendDist/2.0; dr.y = 0.0; dr.z = 0.0;
        printf("placing CMass of dendrimer 1 to -r/2.0,0,0=%lf,%lf,%lf",dr.x,dr.y,dr.z);
        DendPlaceCMass(pd1,dr);
        dr.x = +iniDendDist/2.0; dr.y = 0.0; dr.z = 0.0;
        printf("placing CMass of dendrimer 2 to +r/2.0,0,0=%lf,%lf,%lf",dr.x,dr.y,dr.z);
        DendPlaceCMass(pd2,dr);


/*
 * random vector for dendrimer 1 to dendrimer 2
        RandomVector(&dr,iniDendDist);
        dr.x = pd1->CMass.x + dr.x ;
        dr.y = pd1->CMass.y + dr.y ;
        dr.z = pd1->CMass.z + dr.z ;
*/
        // Set pos equal to npos & calculate centre of Mass
        DendSetpos(pd1);
        DendCalcCMass(pd1);
        
        DendSetpos(pd2);
        DendCalcCMass(pd2);
        
        printf("Position dendrimer 1:"); WriteVecR(pd1->CMass,stdout);
        printf("Position dendrimer 2:"); WriteVecR(pd2->CMass,stdout);


        double dEnergySystem = 0.0;
        dEnergySystem = EnergySystem();
        printf("\nTotal Energy1:%lf\n\n",dEnergySystem);
		FILE *pf;
		pf=fopen("totalpot.dat","w");
		for (double rr=0.5;rr<=1.8;rr+=0.01){
			fprintf(pf,"%13.6f %13.6f\n",rr,Fene(&fnCC,rr)+Morse(&MrsCC,rr));
		}
		fclose(pf);
        MC();
		char filename_lmp[100], filename_coeff[100];
		sprintf(filename_lmp,"in_effdend_r_%08.3f",iniDendDist);
		sprintf(filename_coeff,"coeff_G%dD%d",G,PotType);
		DendWriteinlmp(filename_lmp,filename_coeff,0);
        
        maxRc = MAX3(MrsCC.RCut,MrsCS.RCut,MrsSS.RCut);
        maxRc = 1.1 * maxRc;
        printf("\nMaxRc:%lf\n",maxRc);
        CellListSetup(Box, maxRc);
        CellAlloc();
        CellNewBuild();
        double morse,fene;
        morse  = 0.0;
        fene   = 0.0;
        double result = 0.0;
        struct strMonomer* pmon;
        for (int ndend=0; ndend< numOfDendrimers; ndend++){
            for (int mon=0;mon<dendrimer[ndend].numOfMonomers;mon++){
                pmon    = (dendrimer[ndend].monomer) + mon;
                fene   += FeneEnergy(pmon,&pmon->pos,0);
                morse  += CelldESingleMonMove(pmon,&pmon->pos);
            }
        }
        result = fene/2.0 + morse/2.0;
        printf("Result Cell Energy:%lf\n", result);

    }else{
        printf("Number of dendrimers >2\n");}
    
    return (EXIT_SUCCESS);
}

