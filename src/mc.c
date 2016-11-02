// mc.c
// Added: separate XYZ files 
//        output file with gyration tensor
#include "system.h"
#include "ran250.h"

#define FILE_IDX_LEN(F) (#F)

//#include "preprocessor.h"
void MC(void){
// Open MC file HERE
	int    checkfreq=1;
    double sumEner, sumsqEner;
    
    FILE *pFileDends;
    FILE *pFileEnergy;
	FILE *pFileEquilEnergy;
    FILE *pFileIntraEnergy;
	FILE *pFileGyrTens;
	FILE *pFileCMass;

    pFileDends        = fopen("props.dat","w");
    pFileEnergy       = fopen("energy.dat","w");
	pFileEquilEnergy  = fopen("equilenergy.dat","w");
    pFileIntraEnergy  = fopen("intraenergy.dat","w");
	pFileGyrTens      = fopen("gyrtens.dat","w");
	pFileCMass        = fopen("cmass.dat","w");
//  xyz multiple files declarations
    FILE        *pFileXYZ=NULL;    
	//FILE        *pFileXYZ_PBC;
	int          curfilenum=0,nextfilenum=0; 
	const char  *f_prefix = "ptcls_"; 
	const char  *f_suffix = ".dat"; 
	char         fnameini[120];
	char         fnameiniPBC[120];
	unsigned int numOfFramesFile;
	char        *fname; // allocate space for filenames
	fname    = (char*)malloc(sizeof(char)*(strlen(f_prefix)
			   						     + strlen(FILE_IDX_LEN(numOfXYZfiles))
									     + strlen(f_suffix)));

	if (fname==NULL){
		fprintf(stderr,"no appropriate filename for XYZ specified mc.c\n");
		exit(1);
	}
//  end xyz multiple files
    long unsigned int countTotalSteps=0;
    
    const int adjustcycles=100;
    int i,ncycle;
    // pointer to function to choose between MC_Move 1 or 2
    int (*pMCmove)(int);
    RunningEnergy = EnergySystem();
	printf("MC():: Starting main Monte Carlo Simulation....\n");
    printf("num of displ per cycle:%5d\n",NumofDispPerCycle);
    printf("Running Energy=%lf\n",RunningEnergy);
    EnergyCheck();

    for (i=0;i<3;i++){
        RunningEnergy = EnergySystem3();
        if (i==DECORRELATION){
            NumofCycles = NumofDecorrelationCycles;
			if      (NumofCycles==0)           checkfreq = 1;
			else if ((NumofCycles/100)<=10)    checkfreq = 100;
			else if ((NumofCycles/1000)<=10)   checkfreq = 500;
			else                               checkfreq = NumofCycles/10;
			printf("MC Decorrelation Steps:%-ld\nStart...(check freq=%-7d)\n",NumofCycles,checkfreq);

        }else if (i==EQUILIBRATION){
            NumofCycles = NumofEquilibrationCycles;
			if      (NumofCycles==0)           checkfreq = 1;
			else if ((NumofCycles/100)<=10)    checkfreq = 10;
			else if ((NumofCycles/1000)<=10)   checkfreq = 500;
			else                               checkfreq = NumofCycles/10;
        	if (NumofCycles!=0) printf("MC Equilibration Steps:%-ld\nStart...(check freq=%-7d)\n",NumofCycles,checkfreq);
        }else if (i==PRODUCTION){
			NumofCycles = NumofProductionCycles;
			if      (NumofCycles==0)             checkfreq = 1;
			else if ((NumofCycles/100)<=10)      checkfreq = 100;
			else if ((NumofCycles/1000)<=10)     checkfreq = 5000;
			else if ((NumofCycles/10000)<=10)    checkfreq = 10000;
			else if ((NumofCycles/100000)<=10)   checkfreq = 10000;
			else                                 checkfreq = NumofCycles/100+1;
			if (NumofCycles!=0) printf("MC Production Steps:%-ld\nStart...(check freq=%-7d)\n",NumofCycles,checkfreq);
		}
		EnergyCheck();
		fflush(stdout);
        // set number of attempts and accepted moves to zero
        // for EQUILIBRATION & PRODUCTION
        NumofAttempts      = 0;
        NumofAcceptedMoves = 0;
        //Monte Carlo Cycles

        sumEner    = 0.0;
        sumsqEner  = 0.0;
#if FIXED_CMASS==1
        pMCmove = &McMove_two;
#else
        pMCmove = &McMove;
#endif
        for (ncycle=0;ncycle<NumofCycles;ncycle++){
			//if ((ncycle%checkfreq)==0) {
			//	EnergyCheck();
			//	fflush(stdout);
			//}
#if USE_CELL_LIST==1
			if (numOfDendrimers == 1){
				if (ncycle%60000==0){
					CellNewBuild();
				}
			}
#endif
            countTotalSteps++;

            NumofAcceptedMoves += pMCmove(NumofDispPerCycle);
            NumofAttempts      += NumofDispPerCycle;
            
            sumEner            += RunningEnergy;
            sumsqEner          += RunningEnergy * RunningEnergy;

            if (((ncycle%checkfreq==0))&&(ncycle!=0)){
                double averEner = sumEner/(ncycle+1);
                double devEner  = sqrt((sumsqEner/(ncycle+1)) - SQR(averEner));
                fprintf(stdout,"Step: %d/%ld AverageEnergy: %lf +/- %lf\n",ncycle,NumofCycles,averEner,devEner);
                fflush(stdout);
            }
            if (ncycle%adjustcycles==0){
                AdjustStep(AcceptanceRatio);
            }
            // PRODUCTION MEASUREMENTS
			if (i==EQUILIBRATION || i== DECORRELATION ) {
				if (checkfreq!=0){
					if ((ncycle%checkfreq)==0){
						MCPropEnergy(pFileEquilEnergy);
						fflush(pFileEquilEnergy);
						EnergyCheck();
					}
				}
			}
            if (i==PRODUCTION){
                if ((SamplingFrequency!=0)&&(ncycle%SamplingFrequency==0)){
                    MCProps(pFileDends, pFileEnergy, pFileIntraEnergy, pFileGyrTens);
					DendWriteCMassXYZ(pFileCMass,-1);
                    fflush(pFileDends);
					fflush(pFileIntraEnergy);
                    fflush(pFileEnergy);
					fflush(pFileGyrTens);
					fflush(pFileCMass);
                }
                if ((numOfXYZfiles>0)&&(FramesFrequency!=0)&&(ncycle%FramesFrequency == 0)){
					// !!!! division operator (/) always rounds to zero
					if (ncycle==0) {
						nextfilenum=0;
						numOfFramesFile = NumofCycles/numOfXYZfiles;
					}
					//printf("cur file:%d\n",curfilenum);
					curfilenum = ncycle/numOfFramesFile;
					if ((curfilenum == nextfilenum) ){
						if (curfilenum > 0){
							fclose(pFileXYZ);
							//fclose(pFileXYZ_PBC);
						}
						sprintf(fname,"%s%d%s",f_prefix,curfilenum,f_suffix);
						pFileXYZ = fopen(fname,"w");
						//sprintf(fname,"%sPBC_%d%s\0",f_prefix,curfilenum,f_suffix);
						//pFileXYZ_PBC = fopen(fname,"w");
						nextfilenum++;
					}
					if (nextfilenum > numOfXYZfiles) {
						fprintf(stderr,"mc.c number of XYZ file(%d) exceeds maximum file number(%d)\n",
								nextfilenum,numOfXYZfiles);
					}	
					// Write to XYZ file
                    			DendWriteXYZCommon(pFileXYZ,0);
					sprintf(fnameini,"iniptcls_n%d_G%dF%dD%d_%d.dat",
									 numOfDendrimers,G,F,PotType,Box.x);
					sprintf(fnameiniPBC,"iniptclsPBC_n%d_G%dF%dD%d_%d.dat",
                                        numOfDendrimers,G,F,PotType,Box.x);
					DendSaveConfig(fnameini,0);
					DendSaveConfig(fnameiniPBC,1);
			        //DendWriteXYZCommon(pFileXYZ_PBC,1);
                    fflush(pFileXYZ);
					//fflush(pFileXYZ_PBC);
                }
                
            }

        }
		if (i==DECORRELATION) printf("End of Decorrelation\n");
		if (i==EQUILIBRATION) printf("End of Equilibration\n");
		
        if (i==PRODUCTION){
            printf("*****MC Summary********\n");
            printf("Maxstep:%5.3f\n",MaxStep);
            printf("Acceptance Ratio:%6.3f\n",
                    (double) NumofAcceptedMoves/NumofAttempts);
            printf("Accepted:%ld,Attempts:%ld\n",
                    NumofAcceptedMoves,NumofAttempts);
        }
        EnergyTotal = EnergySystem2();
        double EnergyTotaltest = EnergySystem3();

        printf("Total Energy System:%30.25f Total Energy test:%30.25f Running Energy:%30.25f\n",
                EnergyTotal,EnergyTotaltest,RunningEnergy);
        if (numOfDendrimers == 1) {
            printf("energy dendrimer[0]=%lf , %lf\n",EnergyDend(&dendrimer[0]),IntraDend(&dendrimer[0]));
        }
    }
    EnergyCheck();
    if (numOfXYZfiles!=0){fclose(pFileXYZ);/*fclose(pFileXYZ_PBC);*/}
    fclose(pFileDends);
	fclose(pFileEquilEnergy);
    fclose(pFileEnergy);
    fclose(pFileIntraEnergy);
    fclose(pFileGyrTens);
	fclose(pFileCMass);
}
void MCanneal(double tmax,double tlow, double tstep){ 
	int    isfirstrun;
	int    checkfreq=1;
    double sumEner, sumsqEner;
    
    FILE *pFileDends;
    FILE *pFileEnergy;
	FILE *pFileEquilEnergy;
    FILE *pFileIntraEnergy;
	FILE *pFileGyrTens;
    FILE *pFileXYZ;    

	char         fnameini[120];
	char         fnameiniPBC[120];

	char f_suffix[50];
	char temp_fname[120];

    long unsigned int countTotalSteps=0;
    
    const int adjustcycles=100;
    int i,ncycle;
    // pointer to function to choose between MC_Move 1 or 2
    RunningEnergy = EnergySystem3();
	printf("MC():: Starting main Monte Carlo Simulation....\n");
    printf("num of displ per cycle:%5d\n",NumofDispPerCycle);
    printf("Running Energy=%lf\n",RunningEnergy);
    EnergyCheck();
	// if simulation is restarted then this is NOT
	// the first run => use input number of cycles
	if (Restart == 1) isfirstrun=0;
	else              isfirstrun=1;
	for (Temp=tmax;Temp>=tlow-tstep;Temp-=tstep){
		Beta = 1.0/Temp;
		fprintf(stdout,"Temperature: Temp=%6.3f Beta=%6.3f rho=%6.3lf\n",Temp,Beta,Density);
		sprintf(f_suffix,"_rho%5.3f_T%5.3f.dat",Density,Temp);
		sprintf(temp_fname,"%s%s","props",f_suffix);
		pFileDends        = fopen(temp_fname,"w");
		sprintf(temp_fname,"%s%s","energy",f_suffix);
    	pFileEnergy       = fopen(temp_fname,"w");
		sprintf(temp_fname,"%s%s","equilenergy",f_suffix);
		pFileEquilEnergy  = fopen(temp_fname,"w");
		sprintf(temp_fname,"%s%s","intraenergy",f_suffix);
    	pFileIntraEnergy  = fopen(temp_fname,"w");
		sprintf(temp_fname,"%s%s","gyrtens",f_suffix);
		pFileGyrTens      = fopen(temp_fname,"w");
		sprintf(temp_fname,"%s%s","ptcls",f_suffix);
		pFileXYZ          = fopen(temp_fname,"w");
	    for (i=0;i<3;i++){
				RunningEnergy = EnergySystem3();
				if (i==DECORRELATION){
					if (isfirstrun == 1)
					{
						NumofCycles == 10000;
					}else{
						NumofCycles = NumofDecorrelationCycles;
					}
					if      (NumofCycles==0)           checkfreq = 1;
					else if ((NumofCycles/100)<=10)    checkfreq = 10;
					else if ((NumofCycles/1000)<=10)   checkfreq = 500;
					else                               checkfreq = NumofCycles/10;
        		    printf("MC Decorrelation Steps:%-ld\nStart...(check freq=%-7d)\n",NumofCycles,checkfreq);

        		}else if (i==EQUILIBRATION){
					if (isfirstrun == 1){
						NumofCycles == 510000;
					}else{
						NumofCycles = NumofEquilibrationCycles;
					}
					if      (NumofCycles==0)           checkfreq = 1;
					else if ((NumofCycles/100)<=10)    checkfreq = 10;
					else if ((NumofCycles/1000)<=10)   checkfreq = 500;
					else                               checkfreq = NumofCycles/10;
        		    if (NumofCycles!=0) printf("MC Equilibration Steps:%-ld\nStart...(check freq=%-7d)\n",NumofCycles,checkfreq);
        		}else if (i==PRODUCTION){
        		    NumofCycles = NumofProductionCycles;
					if      (NumofCycles==0)             checkfreq = 1;
					else if ((NumofCycles/100)<=10)      checkfreq = 10;
					else if ((NumofCycles/1000)<=10)     checkfreq = 500;
					else if ((NumofCycles/10000)<=10)    checkfreq = 1000;
					else if ((NumofCycles/100000)<=10)   checkfreq = 10000;
					else                                 checkfreq = NumofCycles/100+1;
        		    if (NumofCycles!=0) printf("MC Production Steps:%-ld\nStart...(check freq=%-7d)\n",NumofCycles,checkfreq);
        		}
				EnergyCheck();
				fflush(stdout);
				// set number of attempts and accepted moves to zero
				// for EQUILIBRATION & PRODUCTION
				NumofAttempts      = 0;
				NumofAcceptedMoves = 0;
				//Monte Carlo Cycles

				sumEner    = 0.0;
				sumsqEner  = 0.0;
				for (ncycle=0;ncycle<NumofCycles;ncycle++){
#if USE_CELL_LIST==1
					if (numOfDendrimers == 1){
						if (ncycle%50000==0){
							CellNewBuild();
						}
					}
#endif
					countTotalSteps++;

					NumofAcceptedMoves += McMove(NumofDispPerCycle);
					NumofAttempts      += NumofDispPerCycle;
					
					sumEner            += RunningEnergy;
					sumsqEner          += RunningEnergy * RunningEnergy;

					if ((((ncycle%checkfreq)==0))&&(ncycle!=0)){
						double averEner = sumEner/(ncycle+1);
						double devEner  = sqrt((sumsqEner/(ncycle+1)) - SQR(averEner));
						fprintf(stdout,"Step: %d/%ld AverageEnergy: %lf +/- %lf\n",ncycle,NumofCycles,averEner,devEner);
						fflush(stdout);
					}
					if (ncycle%adjustcycles==0){
						AdjustStep(AcceptanceRatio);
					}
					// PRODUCTION MEASUREMENTS
					if (i==EQUILIBRATION || i== DECORRELATION ) {
						if (checkfreq!=0){
							if ((ncycle%checkfreq)==0){
								MCPropEnergy(pFileEquilEnergy);
								fflush(pFileEquilEnergy);
							}
						}
					}
					if (i==PRODUCTION){
						if ((SamplingFrequency!=0)&&(ncycle%SamplingFrequency==0)){
							MCProps(pFileDends, pFileEnergy, pFileIntraEnergy, pFileGyrTens);
							fflush(pFileDends);
							fflush(pFileIntraEnergy);
							fflush(pFileEnergy);
							fflush(pFileGyrTens);
						}
						if ( (FramesFrequency!=0)&&(ncycle%FramesFrequency == 0) ){
							// !!!! division operator (/) always rounds to zero
							// Write to XYZ file
							DendWriteXYZCommon(pFileXYZ,0);
							sprintf(fnameini,"iniptcls_n%d_G%dF%dD%d_T%5.3f_%d.dat",
											 numOfDendrimers,G,F,PotType,Temp,Box.x);
							sprintf(fnameiniPBC,"iniptclsPBC_n%d_G%dF%dD%d_T%5.3f_%d.dat",
												numOfDendrimers,G,F,PotType,Temp,Box.x);
							DendSaveConfig(fnameini,0);
							DendSaveConfig(fnameiniPBC,1);
							fflush(pFileXYZ);
						}
						
					}

				}
				if (i==DECORRELATION) printf("End of Decorrelation\n");
				if (i==EQUILIBRATION) printf("End of Equilibration\n");
				
				if (i==PRODUCTION){
					printf("*****MC Summary********\n");
					printf("Maxstep:%5.3f\n",MaxStep);
					printf("Acceptance Ratio:%6.3f\n",
							(double) NumofAcceptedMoves/NumofAttempts);
					printf("Accepted:%ld,Attempts:%ld\n",
							NumofAcceptedMoves,NumofAttempts);
				}
				EnergyTotal = EnergySystem2();
				double EnergyTotaltest = EnergySystem3();

				printf("Total Energy System:%30.25f Total Energy test:%30.25f Running Energy:%30.25f\n",
						EnergyTotal,EnergyTotaltest,RunningEnergy);
				if (numOfDendrimers == 1) {
					printf("energy dendrimer[0]=%lf , %lf\n",EnergyDend(&dendrimer[0]),IntraDend(&dendrimer[0]));
				}
	    }
		EnergyCheck();
    	fclose(pFileXYZ);
    	fclose(pFileDends);
		fclose(pFileEquilEnergy);
    	fclose(pFileEnergy);
    	fclose(pFileIntraEnergy);
    	fclose(pFileGyrTens);
		isfirstrun = 0; // this is not the first run use the number of cycles from
					  // input file
	}

}
void MCequil(int numiter, int numofcycles){
// simple monte carlo routine
// running for numiter times for numofsteps in 
// each iter. 
// Open MC file HERE
    double sumEner, sumsqEner;
    
//  end xyz multiple files
    long unsigned int countTotalSteps=0;
    
    const int adjustcycles=10;
    int i,ncycle;
    // pointer to function to choose between MC_Move 1 or 2
    Beta=1.0/Temp;

	printf("RUNNING EQUILIBRATION MC\n");
	printf("Test Box size:%4d %4d %4d\n",Box.x,Box.y,Box.z);
    RunningEnergy = EnergySystem3();
    printf("num of displ per cycle:%5d\n",NumofDispPerCycle);
    printf("Running Energy=%lf\n",RunningEnergy);
	fflush(stdout);     
    EnergyCheck();
    for (i=0;i<numiter;i++){
		printf("  equil:%d/%d\n",i+1,numiter);
        RunningEnergy = EnergySystem3();
		NumofCycles = numofcycles;
        // set number of attempts and accepted moves to zero
        // for EQUILIBRATION & PRODUCTION
        NumofAttempts      = 0;
        NumofAcceptedMoves = 0;
        //Monte Carlo Cycles

        sumEner    = 0.0;
        sumsqEner  = 0.0;
        for (ncycle=0;ncycle<NumofCycles;ncycle++){
            countTotalSteps++;

            NumofAcceptedMoves += McMove(NumofDispPerCycle);
            NumofAttempts      += NumofDispPerCycle;
            
            sumEner            += RunningEnergy;
            sumsqEner          += RunningEnergy * RunningEnergy;

            if (((ncycle%lround(NumofCycles/10))==0)&(ncycle!=0)){
                double averEner = sumEner/(ncycle+1);
                double devEner  = sqrt((sumsqEner/(ncycle+1)) - SQR(averEner));
                fprintf(stdout,"Step: %d/%d AverageEnergy: %lf +/- %lf\n",ncycle,numofcycles,averEner,devEner);
                //EnergyCheck();
                fflush(stdout);
            }
            if (ncycle%adjustcycles==0){
                AdjustStep(AcceptanceRatio);
            }
            // PRODUCTION MEASUREMENTS
        }
		printf("*****EQUILIBRATION MC Summary********\n");
		printf("Maxstep:%5.3f\n",MaxStep);
		printf("Acceptance Ratio:%6.3f\n",(double) NumofAcceptedMoves/NumofAttempts);
		printf("Accepted:%ld,Attempts:%ld\n",NumofAcceptedMoves,NumofAttempts);
        EnergyTotal = EnergySystem2();
        double EnergyTotaltest = EnergySystem3();

        printf("Total Energy System:%30.25f Total Energy test:%30.25f Running Energy:%30.25f\n",
                EnergyTotal,EnergyTotaltest,RunningEnergy);
    }
    EnergyCheck();
}
int McMove (int NumofDispl){
    int    displ,mtrial,dtrial,accepted;
    Monomer *pm;
    Dendrimer *pd;
#if USE_CELL_LIST!=1
    Dendrimer *pdnext;
#endif
    VecR   *poldPBCpos, newPBCpos, *poldNpos,newNpos;
    
    double energyOld,energyNew,dE;
    
    accepted = 0;
    for(displ=0;displ<NumofDispl;displ++){
        // RANDOM DENDRIMER RANDOM MONOMER
        dtrial = r250n(numOfDendrimers);
        pd = dendrimer+dtrial;
        mtrial = r250n(pd->numOfMonomers);
        pm=(pd->monomer)+mtrial;

        // set old position as the position
        // of the selected monomer
        poldPBCpos = &pm->pos;
        poldNpos   = &pm->npos;

#if USE_CELL_LIST==1
        energyOld = FeneEnergy(pm,poldPBCpos,0) + CelldESingleMonMove(pm,poldPBCpos);
#else
        energyOld = PotMonIntraDend(pm,poldPBCpos,0);
        for (int dend = 0;dend <numOfDendrimers; dend++){
            // test for checking if mon belongs to the dendrimer
            // inside Pot_Mon_NextDend
            pdnext = &dendrimer[dend];
            energyOld += PotMonNextDend(pm,poldPBCpos,pdnext);
        }
        
#endif
        //MonTrialMove(&newPBCpos,poldPBCpos);
        MonTrialMove_All(&newPBCpos,poldPBCpos,&newNpos,poldNpos);

#if USE_CELL_LIST==1
        energyNew = FeneEnergy(pm,&newPBCpos,0) + CelldESingleMonMove(pm,&newPBCpos);
#else

        energyNew   = PotMonIntraDend(pm,&newPBCpos,0);

        for (int dend = 0;dend <numOfDendrimers; dend++){
            pdnext = &dendrimer[dend];
            energyNew += PotMonNextDend(pm,&newPBCpos,pdnext);
        }
#endif


        dE=energyNew-energyOld;

        if (dE<0) {
           // printf("delta epsilon: %21.16f \n",dE);
            poldPBCpos->x = newPBCpos.x;
            poldPBCpos->y = newPBCpos.y;
            poldPBCpos->z = newPBCpos.z;
            poldNpos->x   = newNpos.x;
            poldNpos->y   = newNpos.y;
            poldNpos->z   = newNpos.z;
#if USE_PBC==1
            PBCAll((*poldPBCpos));
#endif
            RunningEnergy+=dE;
            accepted++;
#if USE_CELL_LIST==1
            CellUpdateMon(pm);
#endif
        }

        else if ( (dr250())<(exp(-Beta*dE)) ){
            poldPBCpos->x = newPBCpos.x;
            poldPBCpos->y = newPBCpos.y;
            poldPBCpos->z = newPBCpos.z;
            poldNpos->x   = newNpos.x;
            poldNpos->y   = newNpos.y;
            poldNpos->z   = newNpos.z;
#if USE_PBC==1
            PBCAll((*poldPBCpos));
#endif
            RunningEnergy+=dE;
            accepted++;
#if USE_CELL_LIST==1
            CellUpdateMon(pm);
#endif
        }
    }

    return (accepted);
}
int McMove_two(int NumofDispl){
    // move function for one Monte Carlo Step
    int       displ,mtrial,dtrial,accepted;
    double    energyOld=0.0,energyNew=0.0,energyIntra=0.0,dE=0.0;

    Monomer   *pm;
    Dendrimer *pd;
    Dendrimer *pdnext;

    accepted = 0;

    for(displ=0;displ<NumofDispl;displ++){
        // RANDOM DENDRIMER - RANDOM MONOMER
        dtrial =  r250n(numOfDendrimers);
        pd     =  dendrimer+dtrial;
        mtrial =  r250n(pd->numOfMonomers);
        pm     = &pd->monomer[mtrial];
        // OLD ENERGY
        // energyIntra: summing all internal energy of all OTHER
        // dendrimers that don't contain monomer mtrial
        // it should not change
        // the calculation can be ommited during Monte Carlo steps
        // but not before and after

        energyOld   = IntraDend(pd);

        for (int dend = 0;dend <numOfDendrimers; dend++){
            if ((pd->id)!=dend){
                pdnext       = dendrimer+dend;
                energyOld   += InterDend(pd,pdnext);
            }
        }
        energyOld += energyIntra;
        // MOVE MONOMER
        // keep fixed centre of mass
        MonTrialMoveFixedCMass(pm);
        energyNew = IntraDend(pd);
        for (int dend = 0;dend <numOfDendrimers; dend++){
            if ((pd->id)!=dend){
                pdnext     =dendrimer+dend;
                energyNew += InterDend(pd,pdnext);
            }
        }
        energyNew += energyIntra;
        dE=energyNew-energyOld;

        if (dE<0) {
            RunningEnergy += dE;
            accepted++;
            DendSetnpos(pd);
        }else if (dr250()<exp(-Beta*dE)){
            RunningEnergy += dE;
            accepted++;
            DendSetnpos(pd);
        }else{
            for (int mon=0; mon<pd->numOfMonomers;mon++){
                pd->monomer[mon].pos.x = pd->monomer[mon].bupos.x;
                pd->monomer[mon].pos.y = pd->monomer[mon].bupos.y;
                pd->monomer[mon].pos.z = pd->monomer[mon].bupos.z;
            }
        }
    }

    return (accepted);
}
void AdjustStep(const double accRatio){
    static long unsigned int att=0,acc=0;
    const double acc_min=0.3,acc_max=0.6;
    double curaccratio, curstep;

    if (NumofAttempts==0||att>=NumofAttempts){
        acc = NumofAcceptedMoves;
        att = NumofAttempts;
    }else{
        curaccratio =  (double)(NumofAcceptedMoves-acc)/
                       ((double)(NumofAttempts-att));
// if current acceptance ratio is 30-60% do nothing
// else set maxstep so that the acceptance ration is 
// equal to 50%
		if (curaccratio<acc_min||curaccratio>acc_max){
        	curstep     =  MaxStep;
        	MaxStep     *= fabs(curaccratio/accRatio);
        	if (MaxStep/curstep>1.5)       MaxStep = curstep*1.5;
        	if (MaxStep/curstep<0.5)       MaxStep = curstep*0.5;
        	if (MaxStep>Box.x/4.0) MaxStep = Box.x/4.0;
		}
        	acc = NumofAcceptedMoves;
        	att = NumofAttempts;
    }

}
void AdjustStep_(const double accRatio){
    static long unsigned int att=0,acc=0;
    double curaccratio, curstep;
    if (NumofAttempts==0||att>=NumofAttempts){
        acc = NumofAcceptedMoves;
        att = NumofAttempts;
    }
    else{
        curaccratio =  (double)(NumofAcceptedMoves-acc)/
                       ((double)(NumofAttempts-att));
        curstep     =  MaxStep;
        MaxStep     *= fabs(curaccratio/accRatio);
        if (MaxStep/curstep>1.5)       MaxStep = curstep*1.5;
        if (MaxStep/curstep<0.5)       MaxStep = curstep*0.5;
        if (MaxStep/curstep>Box.x/4.0) MaxStep = Box.x/4.0;
#if DEBUG_MC==1
        fprintf("Old Max Step: %lf\n",curstep);
        fprintf("New Max Step: %lf\n",MaxStep);
        fprintf("Adjusting acc. ratio to:%lf\n",accRatio);
        fprintf("Attempts:%ld\n",NumofAttempts-att);
        fprintf("Accepted:%ld\n",NumofAcceptedMoves-acc);
#endif

    }

}
void MCPropEnergy(FILE *pOutEnergy){
	static unsigned long int count=1;	
    if (count==1){
        fprintf(pOutEnergy,"# [1:frame] [2: EnSys1] [3:EnSys2] [4:EnSys3] "
							"[5:EnRun][6:MorseTotal] [7:FeneTotal] [8:Maxstep] [9:accratio]\n");
    }
    fprintf(pOutEnergy,"%-6lu",count);
    fprintf(pOutEnergy," %16.7e %16.7e %16.7e %16.7e %16.7e %16.7e %13.9f %6.2f\n",
            EnergySystem(),EnergySystem2(),EnergySystem3(),RunningEnergy,MorseTotal,FeneTotal, MaxStep, (double) NumofAcceptedMoves/NumofAttempts);
	count++;
}
void MCProps(FILE *pOutFile, FILE *pOutEnergy,FILE *pOutIntraEnergy,FILE *pOutGyr){
    static unsigned long int count=1;

	double *gyrtens=NULL;
	const short int gyrsize=3;
	const short int elems = SQR(gyrsize);
    int  dendid,i;
    double intraenergy,dummy;
    double DU; //interdendrimer energy 
    //double tol = 1E-5;
    Dendrimer *pd;

//	AllocMat(gyrtens,gyrsize,gyrsize,double);

    if (count==1){
        fprintf(pOutFile,"#[1:Fr.] [2:D.Id] / [3:Intra] / [4:Total] / "
        				 " [5:CMass  x   y    z] [6:Rg] \n");
        fprintf(pOutEnergy,"# [1:frame] [2: EnSys1] [3:EnSys2] [4:EnSys3] "
							"[5:EnRun][6:interdend][7:MorseTotal] [8:FeneTotal] [9:Maxstep] [10:accratio]\n");
		fprintf(pOutGyr,"#1:frame / 3:CM_x / 4:CM_y / 5:CM_z / "
						" 6:x2  7:xy 8:xz 9:yx 10:y2 11:yz 12:zx 13:yz 14:z2\n");
		fprintf(pOutIntraEnergy,"#1:frame #2:columns of intradendrimer energies\n");
    }
    dummy=EnergySystem2();// store energysystem and substract 
						  // from that value the intradend =>interdend
    DU=dummy;
    fprintf(pOutIntraEnergy,"%-6lu",count);
    for (dendid =0; dendid<numOfDendrimers; dendid++){
        pd = dendrimer+dendid;
        DendCalcCMass(pd);
        DendCalcRg(pd);
		gyrtens = DendGyrTensor(pd);
        intraenergy = IntraDend(pd);
		fprintf(pOutIntraEnergy,"  %13.6e",intraenergy);
		DU-=intraenergy;
        fprintf(pOutFile,"%-6lu",count);
        fprintf(pOutFile,"%3d   ",pd->id);
        fprintf(pOutFile,"%13.6e ",intraenergy);
        fprintf(pOutFile,"%13.6e ",RunningEnergy);
        fprintf(pOutFile, "%16.7e %16.7e %16.7e",
                pd->CMass.x, pd->CMass.y, pd->CMass.z);
        fprintf(pOutFile,"%10.3f",pd->Rg.val);
        fprintf(pOutFile,"\n");
		// output gyration file
		fprintf(pOutGyr,"%-6lu %3d %16.7e %16.7e %16.7e "
					   ,count, pd->id, pd->CMass.x, pd->CMass.y, pd->CMass.z);	
		for (i=0;i<elems;fprintf(pOutGyr,"%13.6e ",gyrtens[i]),i++);
		fprintf(pOutGyr,"\n");
    }
    fprintf(pOutIntraEnergy,"\n");
    fprintf(pOutFile,"\n");
	fprintf(pOutGyr,"\n");

    fprintf(pOutEnergy,"%-6lu",count);
    fprintf(pOutEnergy,"   %16.7e %16.7e %16.7e %16.7e %16.7e %16.7e %16.7e",
            EnergySystem(),dummy,EnergySystem3(),RunningEnergy,DU,MorseTotal,FeneTotal);
    fprintf(pOutEnergy," %13.9f %6.2f\n",
            MaxStep, (double) NumofAcceptedMoves/NumofAttempts);
	free(gyrtens);

	count++;
}

