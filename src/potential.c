#include <string.h>
#include <float.h>
#include "system.h"
#include "math.h"

double Fene (FeneParams *fnp,double r){
	if ((r<=fnp->L0+fnp->R)&&(r>=fnp->L0-fnp->R)){
		return -(fnp->K) * SQR(fnp->R)*log(1.0-SQR((r-(fnp->L0))/(fnp->R)));
	}
	else{
		return FLT_MAX;
	}
}
#if LJ==1
static double alj(double r,double rlj){
	double inv_rlj = 1.0/rlj;
	double rinvrlj = r * inv_rlj;
	return 6.0*inv_rlj*(log(rinvrlj)/(rinvrlj-1));
}
#endif
double Morse (MorseParams* mrsp,double r){
#if LJ==1
	double epsLJ=1.0;
	double d=pow(2.0,0.166);
	return epsLJ*(SQR(exp(-alj(r,d)*(r-d))-1.0)-1.0);
#else
    return (mrsp->eps)*(SQR(exp(-(mrsp->a)*(r-(mrsp->d)))-1.0)-1.0);
#endif
}
void PotentialSetup(int dt){
	// 1->7   : correspond to Bianca's models
	// 12     : like D7 but cuttoff and shifted
	// 30->   : customs models mostly for testing purposes
	strcpy(fnCC.type,"CC"); strcpy(fnCS.type,"CS"); strcpy(fnBB.type,"BB");
	strcpy(MrsCC.type,"CC");strcpy(MrsCS.type,"CS");strcpy(MrsSS.type,"SS");
    if (dt==1){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.3; MrsCC.RCut = 2.4;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 2.4;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 2.4;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.7500; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 2.8125; fnBB.R = 0.5625;
    }
    else if (dt==2){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.0; MrsCC.RCut = 2.8;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.0; MrsCS.RCut = 2.8;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 0.0; MrsSS.RCut = 2.8;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.7500; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.875; fnBB.R = 0.3750;
    }
    else if (dt==3){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.45; MrsCC.RCut = 2.6;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.50; MrsCS.RLow = 1.25; MrsCS.RCut = 2.6;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 2.00; MrsSS.RLow = 1.70; MrsSS.RCut = 2.6;
        fnCC.interType  = CC ;fnCC.K = 80.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.75; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 80.0; fnBB.L0 = 2.8125; fnBB.R = 0.5625;
    }
    else if (dt==4){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.45; MrsCC.RCut = 2.6;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.50; MrsCS.RLow = 1.25; MrsCS.RCut = 2.6;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 2.00; MrsSS.RLow = 1.70; MrsSS.RCut = 2.6;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.75; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 2.8125; fnBB.R = 0.5625;
    }
    else if (dt==5){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 4.8;MrsCC.d = 1.00; MrsCC.RLow = 0.25; MrsCC.RCut = 2.7;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.50; MrsCS.RLow = 1.2; MrsCS.RCut = 2.7;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 2.00; MrsSS.RLow = 1.7; MrsSS.RCut = 2.7;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.7500; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 2.8125; fnBB.R = 0.5625;
    }
    else if (dt==6){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.45; MrsCC.RCut = 2.7;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.50; MrsCS.RLow = 1.25; MrsCS.RCut = 2.7;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 2.00; MrsSS.RLow = 1.70; MrsSS.RCut = 2.7;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.7500; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.8750; fnBB.R = 0.375;
    }

    else if (dt==7){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.0; MrsCC.RCut = 2.8;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.0; MrsCS.RCut = 2.8;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 0.0; MrsSS.RCut = 2.8;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 2.8125; fnCS.R = 0.5625;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.8750; fnBB.R = 0.3750;
    }
	else if (dt==11){
        MrsCC.interType = CC ;MrsCC.eps = 0.014;MrsCC.a = 19.2;MrsCC.d = 1.02; MrsCC.RLow = 0.3; MrsCC.RCut = 2.4;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 2.4;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 2.4;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 2.8125; fnCS.R = 0.5625;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.8750; fnBB.R = 0.3750;
	}
	else if (dt==12){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.3; MrsCC.RCut = 1.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 2.4;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 2.4;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 2.8125; fnCS.R = 0.5625;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.8750; fnBB.R = 0.3750;
	}
    else if (dt==30){
		// Similar to D7 but MrsCC.eps = 0.014 instead of 0.714
        MrsCC.interType = CC ;MrsCC.eps = 0.014;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.0; MrsCC.RCut = 2.8;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.0; MrsCS.RCut = 2.8;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 0.0; MrsSS.RCut = 2.8;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 2.8125; fnCS.R = 0.5625;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.8750; fnBB.R = 0.3750;
    }
    else if (dt==40){
		// Similar to D7 but MrsCC.eps = 0.014 instead of 0.714
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 1.80;MrsCC.d = 1.00; MrsCC.RLow = 0.0; MrsCC.RCut = 7.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.01785;MrsCS.a = 6.0;MrsCS.d = 1.75; MrsCS.RLow = 0.0; MrsCS.RCut = 7.0;
        MrsSS.interType = SS ;MrsSS.eps = 0.01785;MrsSS.a = 6.0;MrsSS.d = 2.50; MrsSS.RLow = 0.0; MrsSS.RCut = 7.0;
        fnBB.interType  = BB ;fnBB.K = 60.0; fnBB.L0 = 3.18750; fnBB.R = 0.6375;
        fnCC.interType  = CC ;fnCC.K = 60.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 30.0; fnCS.L0 = 3.5625; fnCS.R = 0.7125;
    }
}


void MorseOutputValues(FILE* outFile){
    fprintf(outFile,"Morse Values\n");
    fprintf(outFile,"CC: rcut=%10.6f phicutoff=%10.6f\n",MrsCC.RCut,Morse(&MrsCC,MrsCC.RCut));
    fprintf(outFile,"CS: rcut=%10.6f phicutoff=%10.6f\n",MrsCS.RCut,Morse(&MrsCS,MrsCS.RCut));
    fprintf(outFile,"SS: rcut=%10.6f phicutoff=%10.6f\n",MrsSS.RCut,Morse(&MrsSS,MrsSS.RCut));

}
void WriteFene(int intertype,double rmin,double rmax,double dr,const char outFileName[]){
    FeneParams *pfn;
    FILE *pOutFile;
    if ((pOutFile=fopen(outFileName,"w"))!=NULL){
        if (intertype==BB) pfn = &fnBB;
        else if (intertype==CC) pfn = &fnCC;
        else if (intertype==CS) pfn = &fnCS;
        else pfn = NULL;
        for (double r=rmin;r<=rmax;r+=dr){
            fprintf(pOutFile,"%lf\t%lf\n",r,Fene(pfn,r));
        }
    }
}
// Simple Potential from array function
void readtoarray(const char infilename[],double x[],double y[],int *numlines){
	int colcounter;
	FILE* pinfile;
	char line[250];	
	char *tmp;
	if ((pinfile=fopen(infilename,"r"))!=NULL){
		while (fgets(line,sizeof line,pinfile)!=NULL){
			//fputs(line,stdout);
			colcounter = 0;
			tmp = strtok(line,"    \t");
			while (tmp !=NULL){
				colcounter++;
				if (colcounter == 1){
					x[*numlines] = strtod(tmp,NULL);
				}
				if (colcounter == 2){
					y[*numlines] = strtod(tmp,NULL);
				}
				//if (colcounter == 3){
				//	y[*numlines] += strtod(tmp,NULL);
				//}
				tmp = strtok(NULL,"\t   ");
			}
			*numlines = *numlines +1;
		}	
	}
}
double readfromarray(double r,double x[],double y[],const int numlines){
    double weight,rbin;
    double mlut,mlut1,mlut2,t1,t2;
	double rmin=x[0];
	double rmax=x[numlines-1];
	double dr=x[2]-x[1];
    int bin;
    double dummy=0.0;
    if ((r<rmax) && (r >= rmin)){
        rbin = (r-rmin)/dr;
        bin  = (int) floor(rbin);
        if (bin <0) {
            printf("readfromarray()::negative bin number\n");
            bin=0;
        }
        weight = rbin - bin;  // rbin is always > rbin (floor)
                              // computer simulations of liquids p.144
        mlut  = y[bin];            // VK
        mlut1 = y[bin+1];          // VK1
        mlut2 = y[bin+2];          // VK2

        t1    = mlut  + (mlut1-mlut) * weight;       //T1
        t2    = mlut1 + (mlut2-mlut1)* (weight-1.0); //T2
        
        dummy = t1 +(t2-t1)*weight*0.5;
    }
    // if r is larger than the cutoff distance of the potential return
    // zero value
    else if (r>rmax)
    {
        dummy = 0.0;
    }

    // if r is less than the minimum distance set potential value
    // to the shortest distance value
    else if (r < rmin){
        dummy = y[0];
    }
    return dummy;
}
//*******************************************//
//*************LOOKUP TABLES*****************//
//*******************************************//
//          FENE Look-up tables              //
//-------------------------------------------//
void LUTMorseAllocate (){
    int interType;
//    type 0: C-C
//    type 1: C-S
//    type 2: S-S
//    0: CC, 1: CS, 2: SS
// allocate morse Table
    if (MorseLUT == NULL){

        MorseLUT        = (double**) malloc ((N_TYPE_MORSE)*sizeof(double*));
        MLUTRsqCut      = (double*)  malloc ((N_TYPE_MORSE)*sizeof(double));
        MLUTRsqMin      = (double*)  malloc ((N_TYPE_MORSE)*sizeof(double));
        MLUTDeltaRsq    = (double*)  malloc ((N_TYPE_MORSE)*sizeof(double));
        MLUTInvDeltaRsq = (double*)  malloc ((N_TYPE_MORSE)*sizeof(double));
        MLUTPhiCutoff   = (double*)  malloc ((N_TYPE_MORSE)*sizeof(double));
        
        for (interType=0;interType<N_TYPE_MORSE;interType++){
            MorseLUT[interType] = (double*) malloc(LUT_MORSE_SIZE*sizeof(double));
        }
    }
}


void LUTMorseBuild (int shift){
//  shift: if 1 then shift potential    
    int interType,bin;
    MorseParams *pmorse;
    for (interType=0;interType<N_TYPE_MORSE;interType++){
    // define parameters according to interaction type
        if      (interType == CC) pmorse=&MrsCC;
        else if (interType == CS) pmorse=&MrsCS;
        else if (interType == SS) pmorse=&MrsSS;
        else{
            printf("LUTMorseBuild():: error "
                    "unknown interaction type %d",interType);
            break;
        }
        MLUTRsqMin[interType]      = SQR(pmorse->RLow);
        MLUTDeltaRsq[interType]    = (SQR(pmorse->RCut) -MLUTRsqMin[interType])/LUT_MORSE_SIZE;
        MLUTInvDeltaRsq[interType] = 1.0/MLUTDeltaRsq[interType];
        MLUTRsqCut[interType]      = SQR(pmorse->RCut);
        MLUTPhiCutoff[interType]   = Morse(pmorse,pmorse->RCut);
        // shift potential 
        double rsq;
        for (bin=0;bin<LUT_MORSE_SIZE;bin++){
            rsq = MLUTRsqMin[interType] + bin*MLUTDeltaRsq[interType];
            MorseLUT[interType][bin] = Morse(pmorse,sqrt(rsq));
			if (shift ==1) MorseLUT[interType][bin]-=MLUTPhiCutoff[interType];
            if (MorseLUT[interType][bin]>DBL_MAX)
                printf("LUTMorseBuild():: MorseLUT[%-1d][%-5d]=%lf", interType,bin,MorseLUT[interType][bin]);
        }
    }
}


double LUTMorseValue (int interType,double rsq){
    
//  Input: interType->interaction type [CC,CS,SS]
//         rsq-> squared distance between particles
//  Output: adds Morse potential value to e
//  if there's no function it adds zero (??????)
//  if rsq < MLUTRsqMin => return value of potential for minimum distance

    double weight,rbin;
    double mlut,mlut1,mlut2,t1,t2;
    int bin;
    double dummy=0.0;
    if ((rsq<MLUTRsqCut[interType]) && (rsq > MLUTRsqMin[interType])){
        rbin = (rsq-MLUTRsqMin[interType])*MLUTInvDeltaRsq[interType];
        bin  = (int) floor(rbin);
        if (bin <0) {
            printf("LUTValue()::negative bin number\n");
            bin=0;
        }
        weight = rbin - bin;// rbin is always > rbin (floor)
        // computer simulations of liquids p.144
        mlut  = MorseLUT[interType][bin];            // VK
        mlut1 = MorseLUT[interType][bin+1];          // VK1
        mlut2 = MorseLUT[interType][bin+2];          // VK2

        t1    = mlut  + (mlut1-mlut) * weight;       //T1
        t2    = mlut1 + (mlut2-mlut1)* (weight-1.0); //T2
        
        dummy = t1 +(t2-t1)*weight*0.5;
/*
        dummy = weight*MorseLUT[interType][bin+1]+
                (1-weight)*MorseLUT[interType][bin];
*/
    }
    // if r is larger than the cutoff distance of the potential return
    // zero value
    else if (rsq>MLUTRsqCut[interType])
    {
        dummy = 0.0;
    }

    // if rsq is less than the minimum distance set potential value
    // to the shortest distance value
    else if (rsq < MLUTRsqMin[interType]){
        dummy = MorseLUT[interType][0];
    }
/*
    printf("LUTValue()::%lf\n",dummy);
*/
    return dummy;
}
//*******************************************//
//          FENE Look-up tables              //
//-------------------------------------------//
void LUTFeneAllocate (){
    int interType;
//    type 0: C-C
//    type 1: C-S
//    type 2: B-B
//    0: CC, 1: CS, 2: BB
    if (FeneLUT == NULL){

        FeneLUT         = (double**) malloc ((N_TYPE_FENE)*sizeof(double*));
        FLUTRsqMin      = (double*)  malloc ((N_TYPE_FENE)*sizeof(double));
        FLUTPhiRsqMin   = (double*)  malloc ((N_TYPE_FENE)*sizeof(double));
        FLUTRsqMax      = (double*)  malloc ((N_TYPE_FENE)*sizeof(double));
        FLUTPhiRsqMax   = (double*)  malloc ((N_TYPE_FENE)*sizeof(double));
        FLUTDeltaRsq    = (double*)  malloc ((N_TYPE_FENE)*sizeof(double));
        FLUTInvDeltaRsq = (double*)  malloc ((N_TYPE_FENE)*sizeof(double));


        for (interType=0;interType<N_TYPE_FENE;interType++){
            FeneLUT[interType] = (double*) malloc(LUT_FENE_SIZE*sizeof(double));
        }
    }
}
void LUTFeneBuild (double factor){
    int interType,bin;
    FeneParams *pfene;
    double rsq,rmin,rmax;
    for (interType=0;interType<N_TYPE_FENE;interType++){
    // define parameters according to interaction type
        if      (interType == CC) pfene=&fnCC;
        else if (interType == CS) pfene=&fnCS;
        else if (interType == BB) pfene=&fnBB;
        else{
            printf("LUTFeneBuild():: error "
                    "unknown interaction type %d",interType);
            break;
        }
        rmin   = (pfene->L0)-factor*(pfene->R);
        rmax   = (pfene->L0)+factor*(pfene->R);
        FLUTRsqMin[interType]      =  SQR(rmin);
        FLUTPhiRsqMin[interType]   =  Fene(pfene,rmin);
        FLUTRsqMax[interType]      =  SQR(rmax);
        FLUTPhiRsqMax[interType]   =  Fene(pfene,rmax);
        FLUTDeltaRsq[interType]    =
                (FLUTRsqMax[interType] -FLUTRsqMin[interType])/LUT_FENE_SIZE;
        FLUTInvDeltaRsq[interType] = 1.0/FLUTDeltaRsq[interType];
        for (bin=0;bin<LUT_FENE_SIZE;bin++){
            rsq = FLUTRsqMin[interType] + bin*FLUTDeltaRsq[interType];
            FeneLUT[interType][bin] = Fene(pfene,sqrt(rsq));
            if (FeneLUT[interType][bin]>DBL_MAX)
                printf("LUTFeneBuild():: FeneLUT[%-1d][%-5d]=%lf",
                        interType,bin,FeneLUT[interType][bin]);
        }
    }
}
double LUTFeneValue (int interType,double rsq){
//  Input: interType->interaction type [CC,CS,BB]
//         rsq-> squared distance between particles
//  Output: adds Fene potential value to e
//  if there's no function it adds zero (??????)
//  if rsq < FLUTRsqMin => return value of potential for minimum distance
    double weight,rbin;
    //double flut,flut1,flut2,t1,t2;
    int bin;
    double dummy=0.0;
    if ((rsq<=FLUTRsqMax[interType]) && (rsq >= FLUTRsqMin[interType])){
        rbin = (rsq-FLUTRsqMin[interType])*FLUTInvDeltaRsq[interType];
        bin  = (int) floor(rbin);
        if (bin <0) {
            printf("LUTValue()::negative bin number\n");
            bin=0;
        }
        weight = rbin - bin;

        //flut  = FeneLUT[interType][bin];            // VK
        //flut1 = FeneLUT[interType][bin+1];          // VK1
        //flut2 = FeneLUT[interType][bin+2];          // VK2

        //t1    = flut  + (flut1-flut) * weight;       //T1
        //t2    = flut1 + (flut2-flut1)* (weight-1.0); //T2
        
        //dummy = t1 +(t2-t1)*weight*0.5;
        dummy = weight*FeneLUT[interType][bin+1]+
                (1-weight)*FeneLUT[interType][bin];
    }
    // if r is larger than the cutoff distance of the potential return
    // lastvalue of table
    else if (rsq >FLUTRsqMax[interType])
    {
//        dummy = FeneLUT[interType][LUT_FENE_SIZE-1];
		dummy = DBL_MAX;
    }
    // if rsq is less than the minimum distance set potential value
    // to the shortest distance table value
    else if (rsq < FLUTRsqMin[interType]){
//        dummy =FeneLUT[interType][0];
		dummy = DBL_MAX;
    }
    return dummy;
}
void LUTMorseOutputTables (const char outFileName[], double rmin, double rmax,double rstep){
    
    FILE  *pOutFile=NULL;
	double pot[N_TYPE_MORSE];
	double potexact[N_TYPE_MORSE];
    MorseParams *pmorse;
	int interType;
    fprintf(stdout,"Testing Morse Lookup Table\n");

    if ((pOutFile=fopen(outFileName,"w"))!=NULL){
		fprintf(pOutFile,"#1:r   2,3:CC   4,5:CS   6,7:SS\n");
		for (double r=rmin;r<rmax;r+=rstep){
			for (interType=0;interType<N_TYPE_MORSE;interType++){
				if      (interType == CC) pmorse=&MrsCC;
        		else if (interType == CS) pmorse=&MrsCS;
        		else if (interType == SS) pmorse=&MrsSS;
        		else{
        		    pmorse = NULL;
        		    printf("LUTMorseOutputTables():: error unknown interaction type %d", interType);
					exit(4);
        		}
				pot[interType]      = LUTMorseValue(interType,SQR(r));
				potexact[interType] = Morse(pmorse,r)-Morse(pmorse,pmorse->RCut);
        	    //fprintf(pOutFile,"%12.6f\t%12.6f\t%12.6f\t",r,e,dif);
        	    //fprintf(pOutFile,"%15.13f\n",dif);
        	}
			fprintf(pOutFile,"%6.3f\t",r);
			for (interType=0;interType<N_TYPE_MORSE;interType++){
				fprintf(pOutFile,"%13.6e  %13.6e\t",pot[interType],potexact[interType]);
			}
			fprintf(pOutFile,"\n");
		}
        
    }
}

void LUTMorseOutputSetup(const char outFileName[]){
    FILE *pOutFile=NULL;
    printf("Testing Morse Look up tables Input\n");
    if ((pOutFile=fopen(outFileName,"w"))!=NULL){
        for (int interType=0; interType<N_TYPE_MORSE;interType++){
            fprintf(pOutFile,"****LUT table properties for type:%-1d*****\n",interType);
            fprintf(pOutFile,"    MLUTDeltaRsq[%1d]%12.8f\n",interType,MLUTDeltaRsq[interType]);
            fprintf(pOutFile,"    MLUTInvDeltaRsq[%1d]%12.8f\n",interType,MLUTInvDeltaRsq[interType]);
            fprintf(pOutFile,"    MLUTPhiCutoff[%1d]%12.8f\n",interType,MLUTPhiCutoff[interType]);
            fprintf(pOutFile,"    MLUTRsqCut[%1d]%12.8f\n",interType,MLUTRsqCut[interType]);
            fprintf(pOutFile,"    MLUTRsqMin[%1d]%12.8f\n",interType,MLUTRsqMin[interType]);
        }
    }
    fclose(pOutFile);
}
void LUTFeneOutputTables (const char outFileName[], double rmin, double rmax,double rstep){

    FILE  *pOutFile=NULL;
	double pot[N_TYPE_FENE];
	double potexact[N_TYPE_FENE];
	int    interType = -1;
    FeneParams *pfene;
    if ((pOutFile=fopen(outFileName,"w"))!=NULL){
		fprintf(pOutFile,"#1:r   2,3:CC   4,5:CS   6,7:BB\n");
		printf("Testing Fene Lookup Table for inter type");
        for (double r=rmin;r<rmax;r+=rstep){
			for (interType=0;interType<N_TYPE_FENE;interType++){
				if      (interType == CC) pfene=&fnCC;
        		else if (interType == CS) pfene=&fnCS;
        		else if (interType == BB) pfene=&fnBB;
        		else{
        		    pfene = NULL;
        		    printf("LUTMorseOutputTables():: error unknown interaction type %d", interType);
        		}
				pot[interType]      = LUTFeneValue(interType,SQR(r));
				potexact[interType] = Fene(pfene,r);
			}
			fprintf(pOutFile,"%6.3f\t",r);
			for (interType=0;interType<N_TYPE_FENE;interType++){
				fprintf(pOutFile,"%13.6e  %13.6e\t",pot[interType],potexact[interType]);
			}
			fprintf(pOutFile,"\n");
        }

    }
}
void LUTFeneOutputSetup(const char outFileName[]){
    FILE *pOutFile=NULL;
    printf("Testing Fene Look up tables Input\n");
    if ((pOutFile=fopen(outFileName,"w"))!=NULL){
        for (int interType=0; interType<N_TYPE_MORSE;interType++){
            fprintf(pOutFile,"****LUT table properties for type:%-1d*****\n",interType);
            fprintf(pOutFile,"    FLUTDeltaRsq[%1d]%12.8f\n",interType,FLUTDeltaRsq[interType]);
            fprintf(pOutFile,"    FLUTInvDeltaRsq[%1d]%12.8f\n",interType,FLUTInvDeltaRsq[interType]);
            fprintf(pOutFile,"    FLUTRsqMin[%1d]%12.8f\n",interType,FLUTRsqMin[interType]);
            fprintf(pOutFile,"    FLUTPhiRsqMin[%1d]%12.8f\n",interType,FLUTPhiRsqMin[interType]);
            fprintf(pOutFile,"    FLUTRsqMax[%1d]%12.8f\n",interType,FLUTRsqMax[interType]);
            fprintf(pOutFile,"    FLUTPhiRsqMax[%1d]%12.8f\n",interType,FLUTPhiRsqMax[interType]);
        }
    }
    fclose(pOutFile);
}
