/* 
 * File:   umbrella.c
 *
 * Created on August 10, 2012, 3:18 PM
 */
#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include <math.h>
#include <time.h>
#include "ran250.h"

typedef unsigned long long int ULONG;

int     HistSize;
double *Hist;
double *HistBackUp;
double *HistRange;
int     lower,mid,upper;
char    outfilename[100];
int     numOfTotalParams=16;

void ReadInput(const char[]);
void Decorrelate(int, int);
void Write(const char outFileName[],int stepfiles,ULONG cycles);

int    DendType;
int    OutputFrequency, CheckFrequency;
double RightBoundary, LeftBoundary, MaximumDisplacement;
double w; // w=f[U(o)]/F[U(n)]  multiplitcation factor for acceptance rule
int    accepted = 0;
double radius;

//static double raninterval(double,double);
void copycoords(Dendrimer*,VecR*,int);
void find_cmass(VecR*,VecR*,int);
void Decorrelate(int ,int);
int Mcmv(int);

int main(int argc, char** argv) {
	time_t time_start,time_cur;
	double time_diff;
	time(&time_start);
    ReadInput(argv[1]);
#ifdef BIAS
	if (argc <3 ){
		printf("Biased potential used, wrong number of input parameters:%d\n",argc -1 );
        exit(EXIT_FAILURE);
	}
	int Numlines=0;
	int Maxlines=6000;
	double *xcol,*ycol;	
	double wn,wo; // rewighting factors
	xcol = (double*)malloc(Maxlines*sizeof(double));
	ycol = (double*)malloc(Maxlines*sizeof(double));

	readtoarray(argv[2],xcol,ycol,&Numlines);	
	
	//FILE *testinpot=fopen("testpot.dat","w");
	//for (double xx=0.0;xx<20.3;xx+=0.01){
	//	fprintf(testinpot,"%lf\t%lf\n",xx,readfromarray(xx,xcol,ycol,Numlines));
	//}
	//fclose(testinpot);
#endif

    if ((PotType>0) && (PotType<=40)){
        PotentialSetup(PotType);
    }else{
        fprintf(stderr,"unknown interaction type\n");
        exit(EXIT_FAILURE);
    }
	FILE *pFileDends       = NULL;
	FILE *pFileEnergy      = NULL;
	FILE *pFileIntraEnergy = NULL;
	FILE *pFileGyrTens     = NULL;
    FILE *pFileXYZ         = NULL;    
	//FILE *pFileAccCom      ;
	if (SamplingFrequency!=0){
		pFileDends        = fopen("props.dat","w");
    	pFileEnergy       = fopen("energy.dat","w");
    	pFileIntraEnergy  = fopen("intraenergy.dat","w");
		pFileGyrTens      = fopen("gyrtens.dat","w");
		pFileXYZ          = fopen("ptcls.dat","w");
	}
	//pFileAccCom = fopen("com_dist.out","w");
    Initialize();
	OutputInputParams();
//  histogram initialization
    Hist=(double*)malloc(HistSize*sizeof(double));
	HistBackUp=(double*)malloc(HistSize*sizeof(double)); 
	HistRange=(double*)malloc((HistSize+1)*sizeof(double));
	for (int h = 0; h <= HistSize; h++)
	{
		double f1 = ((double) (HistSize-h) / (double) HistSize);
		double f2 = ((double) h / (double) HistSize);
		HistRange[h] = f1 * LeftBoundary +  f2 * RightBoundary;
	}
    for (int i=0;i<HistSize;i++) {
		Hist[i]=0.0; 
		HistBackUp[i]=0.0;
	}
//  end histogram initialization
	Decorrelate(0,50000);
	Decorrelate(1,50000);
	DendSetnpos(dendrimer);
	DendSetnpos((dendrimer+1));	
//  Decorrelation step uses only monomer->pos
    VecR vecr0, ranvecr,dr;
	double radius,rr,rr_n;
	
// Give molecules a dendrimer like shape
// place dendrimers inside windows
	vecr0.x = 0.0;
	vecr0.y = 0.0;
	vecr0.z = 0.0;
// move both dendrimers to (0,0,0)
	DendCalcCMass(dendrimer);
	DendCalcCMass(dendrimer+1);
	DendPlaceCMass(dendrimer,vecr0);
	DendPlaceCMass(dendrimer+1,vecr0);
	DendSetpos(dendrimer);
	DendSetpos(dendrimer+1);

	radius = (LeftBoundary+RightBoundary)/2.0;
	RandomVector(&ranvecr,radius);
// move 2nd dendrimer to the middle of the window and
// calculate initial distance
	DendPlaceCMass(dendrimer+1,ranvecr);
	DendSetpos(dendrimer+1);

	printf("initial position:r=%lf\n",VLen(ranvecr));
	WriteVecR(ranvecr,stdout);
	printf("\n");

	VSub(dr,dendrimer->CMass,(dendrimer+1)->CMass);
	rr = VLen(dr);
	if (rr<LeftBoundary||rr>RightBoundary){
		printf("Error... initial molecule distance outside of window\n");
		printf("program exits now\n");
		exit(0);
	}
#if USE_CELL_LIST==1
    maxRc = MAX3(MrsCC.RCut,MrsCS.RCut,MrsSS.RCut);
	CellListSizeMaxrad = 1.6;
    maxRc = CellListSizeMaxrad * maxRc;
    CellListSetup(Box, maxRc);
	CellAlloc();
	CellNewBuild();
#else
	maxRc = -2.0;
#endif
//  MONTE CARLO SIMULATION
    int        mtrial,dtrial;
    Monomer   *pm;
    Dendrimer *pd;
    VecR      *poldPBCpos, newPBCpos, *poldNpos,newNpos;
	ULONG      ncycle;
//BIAS variables
    
    double energyOld,energyNew,dE;
//  temporary positions
	VecR *pos0 = (VecR*)malloc(numOfMonomersInDend*sizeof(VecR));
	VecR *pos1 = (VecR*)malloc(numOfMonomersInDend*sizeof(VecR));	
	VecR   com0 ,com1;
	int    ri;
    

	copycoords(dendrimer  ,pos0,numOfMonomersInDend);		
	copycoords(dendrimer+1,pos1,numOfMonomersInDend);
	find_cmass(pos0,&com0,numOfMonomersInDend);
	find_cmass(pos1,&com1,numOfMonomersInDend);
	VSub(dr,com0,com1);
	rr = VLen(dr);
	if (rr<LeftBoundary||rr>RightBoundary){
		printf("Error... molecule distance outside of window:%lf,[%lf,%lf]\n",rr,LeftBoundary,RightBoundary);
		exit(-1);
	}

	for (int i=0;i<3;i++){
		RunningEnergy  = EnergySystem2();
		EnergyCheck();
		NumofAttempts      = 0;
		NumofAcceptedMoves = 0;
        if (i==DECORRELATION){
            NumofCycles = NumofDecorrelationCycles;
			printf("MC Decorrelation Steps:%-ld\n",NumofCycles);
        }else if (i==EQUILIBRATION){
            NumofCycles = NumofEquilibrationCycles;
			printf("MC Equilibration Steps:%-ld\n",NumofCycles);
        }else if (i==PRODUCTION){
			NumofCycles = NumofProductionCycles;
			printf("MC Production Steps:%-ld\n",NumofCycles);
		}
		fflush(stdout);
		for (ncycle =0; ncycle<NumofCycles;ncycle++){
			if (ncycle%1000==0){
				AdjustStep(0.5);
			}
			for(int displ=0;displ<NumofDispPerCycle;displ++){
	    	    // RANDOM DENDRIMER RANDOM MONOMER
				NumofAttempts++;
	    	    dtrial     = r250n(numOfDendrimers);
				pd         = dendrimer+dtrial;
	    	    mtrial     = r250n(pd->numOfMonomers);
	    	    pm         = (pd->monomer)+mtrial;
        		poldPBCpos = &pm->pos;
        		poldNpos   = &pm->npos;
				// copy coordinates before attempting a move	
#if 	USE_CELL_LIST==1
        		energyOld = FeneEnergy(pm,poldPBCpos,0) + CelldESingleMonMove(pm,poldPBCpos);
#else
        		energyOld = PotMonIntraDend(pm,poldPBCpos,0);
        		for (int dend = 0;dend <numOfDendrimers; dend++){
        		    // test for checking if mon belongs to the dendrimer
        		    // inside Pot_Mon_NextDend
        		    pdnext = dendrimer+dend;
        		    energyOld += PotMonNextDend(pm,poldPBCpos,pdnext);
        		}
#endif
				copycoords(dendrimer  ,pos0,numOfMonomersInDend);		
				copycoords(dendrimer+1,pos1,numOfMonomersInDend);
				find_cmass(pos0,&com0,numOfMonomersInDend);
				find_cmass(pos1,&com1,numOfMonomersInDend);
				VSub(dr,com0,com1);
				rr = VLen(dr);
				
        		MonTrialMove_All(&newPBCpos,poldPBCpos,&newNpos,poldNpos);

				if (dtrial == 0){
					pos0[mtrial].x = newNpos.x;
					pos0[mtrial].y = newNpos.y;
					pos0[mtrial].z = newNpos.z;  
				}else if (dtrial == 1){
					pos1[mtrial].x = newNpos.x;
					pos1[mtrial].y = newNpos.y;
					pos1[mtrial].z = newNpos.z;  
				}
				find_cmass(pos0,&com0,numOfMonomersInDend);
				find_cmass(pos1,&com1,numOfMonomersInDend);
				VSub(dr,com0,com1);
				rr_n = VLen(dr);
	
				if ((rr_n>=LeftBoundary)&&(rr_n<=RightBoundary)){
					//printf("new window distance: %lf\n",rr_n);
#if USE_CELL_LIST==1
			        energyNew = FeneEnergy(pm,&newPBCpos,0) + CelldESingleMonMove(pm,&newPBCpos);
#else
        			energyNew   = PotMonIntraDend(pm,&newPBCpos,0);

        			for (int dend = 0;dend <numOfDendrimers; dend++){
        			    pdnext = dendrimer+dend;
        			    energyNew += PotMonNextDend(pm,&newPBCpos,pdnext);
        			}
#endif
					dE=energyNew-energyOld;
	
#if BIAS==1
					wo=exp(-readfromarray(rr,xcol,ycol,Numlines));
					wn=exp(-readfromarray(rr_n,xcol,ycol,Numlines));
					w=wo/wn;
					//printf("w:%lf wn:%lf wo%lf\n",w,wn,wo);
#else
					w=1.00;
#endif
					if ( (dr250())<(w*exp(-Beta*dE)) ){
						accepted = 1;
					}else{
						accepted = 0;
					}

					if (accepted == 1){
            			poldPBCpos->x  = newPBCpos.x;
            			poldPBCpos->y  = newPBCpos.y;
            			poldPBCpos->z  = newPBCpos.z;
            			poldNpos->x    = newNpos.x;
            			poldNpos->y    = newNpos.y;
            			poldNpos->z    = newNpos.z;
	    	    	    RunningEnergy +=dE;
	    	    	    NumofAcceptedMoves++;
						// Update position
						rr = rr_n;
					    //printf("is accepted: %lf %lf \n",rr,rr_n);
#if 	USE_CELL_LIST==1
						CellUpdateMon(pm);
#endif
					}
					
				}
			}// end of displacements

			if (i==PRODUCTION){
				//printf("rr=%lf, rr_n=%lf\n",rr,rr_n);
				//DendCalcCMass(dendrimer);
				//DendCalcCMass(dendrimer+1);
				//VecR comtest;
				//VSub(comtest,dendrimer->CMass,(dendrimer+1)->CMass);
				//double rrtest=VLen(comtest);
				//printf("calculated distance between com:%lf\n",rrtest);
				//fprintf(pFileAccCom,"%lf\t%lf\n",rr,rrtest);
				//
				//fflush(pFileAccCom);
				upper = HistSize ;
				lower = 0 ;
				
				while (upper - lower > 1)
				{
				mid = (upper + lower) / 2 ;
				
				if (rr >= HistRange[mid])
				{
				  lower = mid ;
				}
				else
				{
				  upper = mid ;
				}
				}
				ri = lower;

				Hist[ri]+=1.0;
			}
			if (CheckFrequency!=0){	
				if (i==PRODUCTION&&(ncycle%CheckFrequency)==0){
					fprintf(stdout,"Cmass , Wrapped CMass:\n");
					VecR tvecr;
					VCopy(tvecr,com0);
					PBCAll(tvecr);
					fprintf(stdout,"Dend[0]:");WriteVecR(com0,stdout);WriteVecR(tvecr,stdout);
					fprintf(stdout,"\n");
					VCopy(tvecr,com1);
					PBCAll(tvecr);
					fprintf(stdout,"Dend[1]:");WriteVecR(com1,stdout);WriteVecR(tvecr,stdout);
					fprintf(stdout,"\n");

					fprintf(stdout,"step:%lld/%ld\n",ncycle,NumofCycles);
					EnergyCheck();
					fprintf(stdout,"Intra 2 dendrimers interactions:%lf\ndistance:%lf\n",IntraDend(dendrimer)+IntraDend((dendrimer+1)),rr);
					fprintf(stdout,"Maximum step:%lf\n",MaxStep);
					time(&time_cur);
					time_diff=difftime(time_cur,time_start);	
					printf("Time elapsed in mins:%lf in hours:%lf\n",time_diff/60,time_diff/3600);
				}
			}
			if (OutputFrequency!=0){
				if (i==PRODUCTION&&(ncycle%OutputFrequency)==0){
					// backup Histogram
					for (int i=0;i<HistSize;i++) {
						HistBackUp[i]=Hist[i];
					}
					Write(outfilename,0,ncycle);	
				}
			}
			if (SamplingFrequency!=0){
				if (i==PRODUCTION&&(ncycle%SamplingFrequency)==0){
					MCProps(pFileDends, pFileEnergy, pFileIntraEnergy, pFileGyrTens);
					DendWriteXYZCommon(pFileXYZ, 0);
                	fflush(pFileDends);
					fflush(pFileIntraEnergy);
                	fflush(pFileEnergy);
					fflush(pFileGyrTens);
					fflush(pFileXYZ);
				}
			}
		}
        if (i==DECORRELATION||i==PRODUCTION){
            printf("*****MC Summary********\n");
            printf("Maxstep:%5.3f\n",MaxStep);
            printf("Acceptance Ratio:%6.3f\n",
                    (double) NumofAcceptedMoves/NumofAttempts);
            printf("Accepted:%ld,Attempts:%ld\n",
                    NumofAcceptedMoves,NumofAttempts);
        }
	}
	if (SamplingFrequency!=0){
		fclose(pFileDends);
    	fclose(pFileEnergy);
    	fclose(pFileIntraEnergy);
    	fclose(pFileGyrTens);
		fclose(pFileXYZ);
	}
	//fclose(pFileAccCom);
	
	time(&time_cur);
	time_diff=difftime(time_cur,time_start);	
	printf("Total time elapsed in mins:%lf in hours:%lf\n",time_diff/60,time_diff/3600);
    return (EXIT_SUCCESS);
}
void copycoords(Dendrimer* pd,VecR *pos,int N){
	if (N!=pd->numOfMonomers) {
		printf("error while copying coords... \nnum of monomers in dend differs from input\n");
		exit(0);
	}
    for (int mon = 0; mon < pd->numOfMonomers; mon++) {
        pos[mon].x = (pd)->monomer[mon].npos.x;
        pos[mon].y = (pd)->monomer[mon].npos.y;
        pos[mon].z = (pd)->monomer[mon].npos.z;
    }
	
}
void find_cmass(VecR *pos,VecR *com,int N){

	com->x = 0.0;
	com->y = 0.0;
	com->z = 0.0;

	for (int i=0;i<N;i++){
		com->x += pos[i].x;
		com->y += pos[i].y;
		com->z += pos[i].z;
	}	
	com->x /= (double)N;
	com->y /= (double)N;
	com->z /= (double)N;
}

void Write(const char outFileName[],int stepfiles,ULONG cycles){
	static int fcounter=0;
    int ri;
	double rr,vol,dummy,dummybak,deltaR;
	FILE *poutfile=NULL,*poutbackfile=NULL; 
	char filename[200];
	if ((stepfiles==1)&&(fcounter<10000)&&(cycles%200000==0)) {
		sprintf(filename,"b%05d.bakh",fcounter);
		poutbackfile = fopen(filename,"w");
	}
		
    if ((poutfile=fopen(outFileName,"w"))!=NULL){
		fprintf(poutfile,"# number of cycles:%lld\n",cycles);
		fprintf(poutfile,"#1:rindex, 2:r 3:rwmin 4:rwmax, 5: hist_value 6: numof cycles 7:shell volume 8: normalized hist value 9: hist_values/ number of cycles 10: -log(hist_value/numcycles)\n");
    	for (ri=0;ri<HistSize;ri++){
			rr     = (HistRange[ri]+HistRange[ri+1])/2.0; // shell radius
			deltaR = HistRange[ri+1]-HistRange[ri]; // shell dr
			vol    = 4*PI*SQR(rr)*deltaR;//calculate shell volume 4*pi*r^2*dr
    	    dummy  = Hist[ri]/vol;
    	    fprintf(poutfile,"%3d\t%13.6E\t%13.6E\t%13.6E\t%13.6E\t%lld\t%13.6E\t%13.6E\t%13.6E\t%13.6E\n",
							ri,rr,HistRange[ri],HistRange[ri+1],Hist[ri],cycles,vol,dummy,dummy/cycles,-log(dummy/cycles));
			if ((stepfiles == 1 )&&(fcounter<10000)&&(cycles%200000==0)){
				dummybak = HistBackUp[ri]/vol; 
				fprintf(poutbackfile,"%3d\t%13.6E\t%13.6E\t%13.6E\t%13.6E\t%lld\t%13.6E\t%13.6E\t%13.6E\t%13.6E\n",
							ri,rr,HistRange[ri],HistRange[ri+1],HistBackUp[ri],cycles,vol,dummybak,dummybak/cycles,-log(dummybak/cycles));
			}
    	}
	}
	fflush(poutfile);
	fclose(poutfile);
	if (stepfiles ==1){
		fflush(poutbackfile);
		fclose(poutbackfile);
	}
	fcounter++;
}


void ReadInput(const char inFileName[]){
    int count=0; // count parameters read
	int Boxl;
    numOfDendrimers = 2;
    seed=0;
    Temp = 1.0;
    MaxStep=0.0;
    numOfIniConfigFiles = 0;
    FILE *pInFile;
    // file format:
    if ((pInFile=fopen(inFileName,"r"))!=NULL){
		printf("Reading Params\n");
		count+=fscanf(pInFile,"%d\n",&Boxl);
        count+=fscanf(pInFile,"%lf\n",&LeftBoundary);
        count+=fscanf(pInFile,"%lf\n",&RightBoundary);
        count+=fscanf(pInFile,"%d\n" ,&HistSize);
		count+=fscanf(pInFile,"%s"   ,outfilename);
        count+=fscanf(pInFile,"%d\n" ,&PotType);
        count+=fscanf(pInFile,"%d\n" ,&G);
        count+=fscanf(pInFile,"%d\n" ,&F);
        count+=fscanf(pInFile,"%ld\n",&NumofDecorrelationCycles);
        count+=fscanf(pInFile,"%ld\n",&NumofEquilibrationCycles);
        count+=fscanf(pInFile,"%ld\n",&NumofProductionCycles);
        count+=fscanf(pInFile,"%u\n" ,&seed);
        count+=fscanf(pInFile,"%lf\n",&MaxStep);
		count+=fscanf(pInFile,"%d\n" ,&OutputFrequency);
		count+=fscanf(pInFile,"%d\n" ,&CheckFrequency);
		count+=fscanf(pInFile,"%d\n" ,&SamplingFrequency);
    }
	else{
		printf("Cannot find input file:%s\n",inFileName);
	}
	Box.x =Boxl; Box.y=Boxl;Box.z = Boxl;
	fprintf(stdout,"Box Size:%d\n",Box.x);
	fprintf(stdout,"Left Boundary:%lf\n",LeftBoundary);
	fprintf(stdout,"Right Boundary:%lf\n",RightBoundary);
	fprintf(stdout,"histogram size:%d\n", HistSize);
	fprintf(stdout,"output filename: %s\n",outfilename);
	fprintf(stdout,"dendrimer type: %d\n", PotType);
	fprintf(stdout,"generation:%d\n", G);
	fprintf(stdout,"functionality:%d\n", F);
	fprintf(stdout,"decorrelation cycles:%ld\n",NumofDecorrelationCycles);
	fprintf(stdout,"equilibration cycles:%ld\n",NumofEquilibrationCycles);
	fprintf(stdout,"production cycles:%ld\n",NumofProductionCycles);
	fprintf(stdout,"seed:%u\n",seed);
	fprintf(stdout,"output every %d steps\n", OutputFrequency);
	fprintf(stdout,"check  every %d steps\n", CheckFrequency);
	fprintf(stdout,"sample  every %d steps\n", SamplingFrequency);
	    if (count!=numOfTotalParams){
        fprintf(stderr,"num of params read:%2d, expected:%2d\n",count,numOfTotalParams);
    }
}
//static double raninterval(double dblmin,double dblmax){
//    // returns a random double number from (dblmin,dblmax)
//    return ((dblmax-dblmin)*dr250() + dblmin );
//}
void Decorrelate(int dendid, int tries){
    for (int i=0;i<tries;i++){
        NumofAttempts      += tries;
        NumofAcceptedMoves += Mcmv(dendid);
    }
}
int Mcmv(int dendid){
    int        displ,mtrial,accepted;
    Monomer   *pm;
    Dendrimer *pd;
    VecR      *poldPos,newPos;
    double     energyOld,energyNew,dE=0.0;
    int        NumofDispl;
    accepted = 0;

    pd = dendrimer+dendid;
    NumofDispl = pd->numOfMonomers;

    for(displ=0;displ<NumofDispl;displ++){
        // RANDOM DENDRIMER RANDOM MONOMER
        mtrial = r250n(pd->numOfMonomers);
        pm     = (pd->monomer)  + mtrial;
        // set old position as the position
        // of the selected monomer
        poldPos=&pm->pos;

        energyOld = PotMonIntraDend(pm,poldPos,0);

        MonTrialMove(&newPos,poldPos);

        energyNew = PotMonIntraDend(pm,&newPos,0);


        dE=energyNew-energyOld;

        if (dE<0) {
            poldPos->x=newPos.x;
            poldPos->y=newPos.y;
            poldPos->z=newPos.z;
            RunningEnergy+=dE;
            accepted++;
        }

        else if ( (dr250())<(exp(-dE)) ){
            poldPos->x=newPos.x;
            poldPos->y=newPos.y;
            poldPos->z=newPos.z;
            RunningEnergy+=dE;
            accepted++;
        }
    }
     
    return (accepted);
}
