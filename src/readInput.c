// readInput.c 
#include "system.h"
void ReadSimParamsAnneal(const char inFileName[]){
    FILE *pInFile;
    int  count=0;
    if ((pInFile=fopen(inFileName,"r"))!=NULL){
        count+=fscanf(pInFile,"%d %d %d\n",&G,&F,&numOfDendrimers);
        count+=fscanf(pInFile,"%u %d %s %d\n",&seed,&numOfIniConfigFiles,fileIniConf,&Restart);
        count+=fscanf(pInFile,"%lf %lf %lf\n",&TempHigh,&TempLow,&TempStep);
        count+=fscanf(pInFile,"%lf\n",&Temp);
        count+=fscanf(pInFile,"%d %d %d\n",&Box.x,&Box.y,&Box.z);
		count+=fscanf(pInFile,"%lf\n",&CellListSizeMaxrad);
        count+=fscanf(pInFile,"%ld %ld %ld %lf\n",&NumofDecorrelationCycles,&NumofEquilibrationCycles,&NumofProductionCycles,&AcceptanceRatio);
        count+=fscanf(pInFile,"%d %d\n",&SamplingFrequency,&FramesFrequency);
        count+=fscanf(pInFile,"%lf\n",&MaxStep);
		count+=fscanf(pInFile,"%d %s\n",&PotType,filePotential);
	} else{
        printf("input file (%s) not found\n",inFileName);
        exit(3);
    }
}


void ReadSimParams(const char inFileName[]){
    const int inputparams = 30;
    FILE *pInFile;
    int  count=0;
    if ((pInFile=fopen(inFileName,"r"))!=NULL){

		count+=fscanf(pInFile,"%d",&TaskID);
        count+=fscanf(pInFile,"%d %d %d\n",&G,&F,&numOfDendrimers);
        count+=fscanf(pInFile,"%lf\n",&iniDendDist);
        count+=fscanf(pInFile,"%u %d %s %d\n",&seed,&numOfIniConfigFiles,fileIniConf,&Restart);
		count+=fscanf(pInFile,"%d\n",&numOfXYZfiles);
        count+=fscanf(pInFile,"%lf %d %d %d\n",&Temp,&Box.x,&Box.y,&Box.z);
		count+=fscanf(pInFile,"%lf\n",&CellListSizeMaxrad);
        count+=fscanf(pInFile,"%ld %ld %ld %lf\n",&NumofDecorrelationCycles,&NumofEquilibrationCycles,&NumofProductionCycles,&AcceptanceRatio);
        count+=fscanf(pInFile,"%d %d\n",&SamplingFrequency,&FramesFrequency);
        count+=fscanf(pInFile,"%lf\n",&MaxStep);
		count+=fscanf(pInFile,"%d %s\n",&PotType,filePotential);
        if (count!=inputparams){
            printf("wrong number of input params...\nexpected:%-4d\tactual:%-4d\n",
                    inputparams,count);
//            exit(3);
        }
       
    }
    else{
        printf("input file (%s) not found\n",inFileName);
        exit(3);
    }
    
}

void ReadPotential(int pottype,const char inFileName[]){
	//const int inputparams = 30;
    FILE *pInFile;
    int  count=0;
	if (pottype==0){
			strcpy(fnCC.type,"CC"); strcpy(fnCS.type,"CS"); strcpy(fnBB.type,"BB");
			strcpy(MrsCC.type,"CC");strcpy(MrsCS.type,"CS");strcpy(MrsSS.type,"SS");
			if ((pInFile=fopen(inFileName,"r"))!=NULL){
				count+=fscanf(pInFile,"%d %lf %lf %lf\n",
							  &fnCC.interType,&fnCC.K,&fnCC.L0,&fnCC.R);
				count+=fscanf(pInFile,"%d %lf %lf %lf\n",
							  &fnCS.interType,&fnCS.K,&fnCS.L0,&fnCS.R);
				count+=fscanf(pInFile,"%d %lf %lf %lf\n",
							  &fnBB.interType,&fnBB.K,&fnBB.L0,&fnBB.R);
#if LJ!=1
				count+=fscanf(pInFile,"%d %lf %lf %lf %lf %lf\n",
							  &MrsCC.interType,&MrsCC.eps,&MrsCC.a,&MrsCC.d,&MrsCC.RLow,&MrsCC.RCut);
				count+=fscanf(pInFile,"%d %lf %lf %lf %lf %lf\n",
							  &MrsCS.interType,&MrsCS.eps,&MrsCS.a,&MrsCS.d,&MrsCS.RLow,&MrsCS.RCut);
				count+=fscanf(pInFile,"%d %lf %lf %lf %lf %lf\n",&MrsSS.interType,
							  &MrsSS.eps,&MrsSS.a,&MrsSS.d,&MrsSS.RLow,&MrsSS.RCut);
#else
				double dummyRCut = pow(2,0.1666);
				MrsCC.RCut = dummyRCut;
				MrsCS.RCut = dummyRCut;
				MrsSS.RCut = dummyRCut;
#endif
			}
	}else{
		PotentialSetup(pottype);	
	}
}

void OutputInputParams(){
        fprintf(stdout,"*******Input Parameters*******\n");
        fprintf(stdout,"------------------------------\n");
        fprintf(stdout,"# of dendrimers:%d\n# of Generations: %-d\n"
                "Functionality:%-d\n",numOfDendrimers,G,F);
        fprintf(stdout,"Temperature:%-6.2f\nBox::  x=%-5.1d y=%-5.1d z=%-5.1d\n",
                Temp,Box.x,Box.y,Box.z);
     	fprintf(stdout,"num of XYZ particles files:%d\n",numOfXYZfiles);
		if (Restart == 1) fprintf(stdout,"Restarting simulation\n");
		else              fprintf(stdout,"No restart\n");
     	fprintf(stdout,"CELL LIST cell size %lf*maxrad\n",CellListSizeMaxrad);
        fprintf(stdout,"MC parameters\n");
        fprintf(stdout,"# of Equilibration cycles=%-ld\n# of Production cycles %-ld\n",
                NumofEquilibrationCycles,NumofProductionCycles);
        fprintf(stdout,"Sampling Frequency=%-4d\n# of Frames VMD: %-4d\n",
                SamplingFrequency,FramesFrequency);
        fprintf(stdout,"Maxstep:%5.3f\n",MaxStep);
		fprintf(stdout,"Potential Type:%d %s\n",PotType,filePotential);	
        fprintf(stdout,"Fene Potential Paramaters\n");
        fprintf(stdout,"type|   K   |   L0   |   R   |\n");
        fprintf(stdout,"-------------------------\n");
        fprintf(stdout,"%3s  %7.3f %7.3f %7.3f\n",fnCC.type,fnCC.K,fnCC.L0,fnCC.R);
        fprintf(stdout,"%3s  %7.3f %7.3f %7.3f\n",fnCS.type,fnCS.K,fnCS.L0,fnCS.R);
        fprintf(stdout,"%3s  %7.3f %7.3f %7.3f\n",fnBB.type,fnBB.K,fnBB.L0,fnBB.R);
        fprintf(stdout,"\n");
        fprintf(stdout,"Morse Potential Paramaters\n");
        fprintf(stdout,"type|  eps  |   a    |   d   | Rlow | RCut |\n");
        fprintf(stdout,"--------------------------------------------\n");
        fprintf(stdout,"%3s %7.3f %7.3f %7.3f %7.3f %7.3f\n",
                MrsCC.type,MrsCC.eps,MrsCC.a,MrsCC.d,MrsCC.RLow,MrsCC.RCut);
        fprintf(stdout,"%3s %7.3f %7.3f %7.3f %7.3f %7.3f\n",
                MrsCS.type,MrsCS.eps,MrsCS.a,MrsCS.d,MrsCS.RLow,MrsCS.RCut);
        fprintf(stdout,"%3s %7.3f %7.3f %7.3f %7.3f %7.3f\n\n",
                MrsSS.type,MrsSS.eps,MrsSS.a,MrsSS.d,MrsSS.RLow,MrsSS.RCut);
        fprintf(stdout,"Initial conf file name :%s\n\n",fileIniConf);
        fprintf(stdout,"num of initial config file:%-3d\n",numOfIniConfigFiles);
        fflush(stdout);
}
