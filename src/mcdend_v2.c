/* 
 * File:   mcdend.c
 * Author: georgiou
 *
 * Created on February 9, 2011, 10:57 AM
 */
#include "system.h"
#include "ran250.h"

int main() {
    int i,j;
    ReadInputData();
    SetParams();
    OpenOutFiles();
    r250_init(seed);
    AllocArrays();
#if USELUT==1
    printf("Using LookupTables\n");
    DefineLookUpTables();
#endif
    Initialize();
    InitConfig();
    DefineBonds();
    EnergySystem();
    if (EnergyTotal>DBL_MAX) {
        printf("Error in Energy of initial configuration");
        exit(0);
    }
    WriteBondVmd(pBonds);
    RunningEnergy = EnergyTotal;
    printf("EnergyTotal before Monte Carlo: %-100.6f\n",EnergyTotal );
    
    for (i=0;i<2;i++){
        if (i==0){
            NumofCycles = NumofEquilibrationCycles;
            if (NumofCycles!=0)
                printf("MC Equilibration Steps:%-10d\nStart...\n",NumofCycles);
        } else{
            NumofCycles = NumofProductionCycles;
            if (NumofCycles!=0) printf("MC Production Steps:%-10d\nStart...\n",NumofCycles);
        }

        NumofAttempts            = 0;
        NumofAcceptedMoves       = 0;
        // Initializing AdjustStep
        AdjustStep(0.5);

        for (j=0;j<NumofCycles;j++){
            if ((i==0)&&(j==100)){
                printf("Initializing Total Energy\n");
                EnergySystem();
                RunningEnergy=EnergyTotal;
                printf("Total Energy=%-50.7f\n",EnergyTotal);
            }
            NumofAttempts      += NumofDispPerCycle;
            // MC Step
            NumofAcceptedMoves += McMove(NumofDispPerCycle);

            if ((j%(NumofCycles/5)) == 0){
                AdjustStep(0.5);
            }
            if (i==1){
                if ((j*NumofDispPerCycle%1000000)==0)
                    printf("MC Step: %d/%d\n",j,NumofCycles);
                if ((j%SamplingFrequency) == 0){
                    CalcProps(pProps,pCMass);
                }
                if(j%(NumofCycles/NumofFrames)==0){
                    WriteXYZ(pMC);
                }
            }
        }
        if (i==1){
            AcceptanceRatio = (double) NumofAcceptedMoves/NumofAttempts;
            printf("Acceptance Ratio:%6.3f\n",AcceptanceRatio);
            printf("Accepted:%10d,Attempts:%10d\n",NumofAcceptedMoves,NumofAttempts);
        }
        EnergySystem();
        printf("Total Energy System:%-lf\n",EnergyTotal);
    }
    printf("Run Energy MC:      %-lf\n",RunningEnergy);
    
    CloseOutFiles();
}

void ReadInputData(){
    // Seed for random number generator
    seed=65456091;
    // Dendrimer properties
    G=6;
    P=1;
    F=3;
    // Maximum displacement for MC step
    // The max. displ. is adjusted with function Adjust()
    MaxStep=0.3;


    // Fene Parmeters
    strcpy(fnCC.type,"CC"); fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
    strcpy(fnCS.type,"CS"); fnCS.K = 30.0; fnCS.L0 = 3.7500; fnCS.R = 0.7500;
    strcpy(fnBB.type,"BB"); fnBB.K = 40.0; fnBB.L0 = 2.8125; fnBB.R = 0.5625;
    // Morse Parameters
    strcpy(MrsCC.type,"CC");MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00;
    strcpy(MrsCS.type,"CS");MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.50;
    strcpy(MrsSS.type,"SS");MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 2.00;

    RCut = 2.70;

    // Look up tables definitions
#if USELUT==1
    TableSize = 200;
    TableRmin = 0.0;
    RCutoffCC =1.6;RCutoffCS =1.5;RCutoffSS =1.96;
#endif
    // MC Params
    NumofEquilibrationCycles = 100000;
    NumofProductionCycles    = 1200000;
    NumofFrames              = 2000;
    
    SamplingFrequency        = 80;

}

void SetParams(void){
    
    double cc;
    int L;

    Nmon = 0;

    for (int ngen=0;ngen<=G;ngen++){
        Nmon +=2*(int) pow((F-1),ngen);
        if (ngen==(G-1)) NCmon = Nmon;
    }
/*
    Nmon = (F-1)*(pow((double)(F-1),(double)(G+1)));
*/
    NSmon = Nmon - NCmon;
    printf("Number of Monomers: %-4d \n",Nmon);
    printf("# core monomers:%-3d\n# shell monomers:%-3d\n",NCmon,NSmon);
    
    // Range and number of bins of the histogram
    dpSizeHist =200;
    // Determine Histogram Range;
    dpRange = ceil(2*G*MAX3(fnBB.L0,fnCC.L0,fnCS.L0)+fnBB.L0);
#if DEBUG == 0
    printf("*****Histogram Properties******\n");
    printf("Range of Histogram:=%-8.3f\n",dpRange);
    printf("Number of Bins:=%-d\n",dpSizeHist);
#endif
    // Box size
    L=LPM*Nmon;
    cc=1.0/NDIM;
    L=LPM*Nmon;
    L=pow(L,cc);
    Box.x=L+1;Box.y=L+1;Box.z=L+1;

    // Set # of displacements / MC Step
    NumofDispPerCycle = Nmon;
    FrameStep         = NumofProductionCycles/NumofFrames;
    // Dp histogram limit
    dpLimit =NumofProductionCycles / SamplingFrequency;
    
}

void   OpenOutFiles(){
    char fnIniConfig[40],fnMC[40],fnEnergy[40],fnBonds[40],fnCMass[40];
    sprintf(fnIniConfig,"iniconfig_G%dF%d.xyz",G,F);
    sprintf(fnMC,"MC_G%dF%d.xyz",G,F);
    sprintf(fnEnergy,"EnRg_G%dF%d.xyz",G,F);
    sprintf(fnBonds,"Bonds_G%dF%d.vmd",G,F);
    sprintf(fnCMass,"CMass_G%dF%d.xyz",G,F);
    pIniConfig = fopen(fnIniConfig,"w");
    pMC        = fopen(fnMC,"w");
    pBonds     = fopen(fnBonds,"w");
    pProps     = fopen(fnEnergy,"w");
    pCMass     = fopen(fnCMass,"w");
}
void   CloseOutFiles(){
    fclose(pIniConfig);
    fclose(pMC);
    fclose(pBonds);
    fclose(pProps);
    fclose(pCMass);
}

void AllocArrays(){
    Mon = (Monomer*) malloc(Nmon*sizeof(Monomer));
    for (int m=0;m<Nmon;m++){
        Mon[m].bond  = (int*) malloc(F*sizeof(int));
    }
    // Allocate DP histogram
    dpHist=(double**) malloc(NTYPE*sizeof(double *));
    for (int n=0;n<NTYPE;n++){
        dpHist[n]=(double*)malloc(dpSizeHist*sizeof(double));
    }
    dpHistGen=(double**) malloc ((G+1)*sizeof(double*));
    for (int ngen=0;ngen<=G;ngen++){
        dpHistGen[ngen]=(double*)malloc(dpSizeHist*sizeof(double));
    }
#if USELUT==1
    // Allocate LUT
    MorseTab = (double**) malloc ((LUT_NINTER)*sizeof(double*));
    for (int interType=0;interType<LUT_NINTER;interType++){
        MorseTab[interType] = (double*) malloc(TableSize*sizeof(double));
    }
    TableDeltaRsq    = (double*) malloc ((LUT_NINTER)*sizeof(double));
    TableInvDeltaRsq = (double*) malloc ((LUT_NINTER)*sizeof(double));
    TableRCutOff     = (double*) malloc ((LUT_NINTER)*sizeof(double));
    TableMorseCutoff = (double*) malloc ((LUT_NINTER)*sizeof(double));
#endif

}
#if USELUT
void   DefineLookUpTables (void){
    // Compute once tables for morse potential
    MorseParams *pMrsPrms;
    int interType,bin;
    double rsq;
    TableRsqMin = SQR(TableRmin);

    for (interType=0;interType<LUT_NINTER;interType++){
        // define parameters according to interaction type
        if (interType == CC) {
            pMrsPrms=&MrsCC;
            TableDeltaRsq[CC]    = (SQR(RCutoffCC)-TableRsqMin)/TableSize;
            TableInvDeltaRsq[CC] = 1.0/TableDeltaRsq[CC];
            TableRCutOff[CC]     = SQR(RCutoffCC);
            TableMorseCutoff[CC] = Morse(pMrsPrms,RCutoffCC);
        }
        else if (interType == CS) {
            pMrsPrms=&MrsCS;
            TableDeltaRsq[CS]    = (SQR(RCutoffCC)-TableRsqMin)/TableSize;
            TableInvDeltaRsq[CS] = 1.0/TableDeltaRsq[CS];
            TableRCutOff[CS]     = SQR(RCutoffCS);
            TableMorseCutoff[CS] = Morse(pMrsPrms,RCutoffCS);
        }
        else if (interType == SS){
            pMrsPrms=&MrsSS;
            TableDeltaRsq[SS]=(SQR(RCutoffCC)-TableRsqMin)/TableSize;
            TableInvDeltaRsq[SS] = 1.0/TableDeltaRsq[SS];
            TableRCutOff[SS]     = SQR(RCutoffSS);
            TableMorseCutoff[SS] = Morse(pMrsPrms,RCutoffSS);
        }
        else{
            printf("DefineLookUpTables():: error "
                    "unknown interaction type %d",interType);
            break;
        }

        for (bin=0;bin<TableSize;bin++){
            rsq = TableRsqMin + bin*TableDeltaRsq[CC];
            MorseTab[interType][bin] = Morse(pMrsPrms,sqrt(rsq))
                    -TableMorseCutoff[interType];
        }
    }

}
#endif

void Initialize (){
    for (int m=0;m<Nmon;m++){
        Mon[m].gen=-1;
        for (int nbond=0;nbond<F;nbond++){
            Mon[m].bond[nbond]=-2;
        }
    }
}
void DefineBonds (){
    // Set properties for each monomer
    // id, generation, bonds,type
    // Central Monomers
    int Ngf,Ngl,mnext,nbond;  // number of first and last monomer of current generation
    Mon[0].id        = 0;
    Mon[0].type      = B;
    Mon[0].gen       = 0;
    Mon[0].bond[F-1] = 1;

    Mon[1].id        = 1;
    Mon[1].type      = B;
    Mon[1].gen       = 0;
    Mon[1].bond[F-1] = 0;

    Ngf=0;
    mnext = 2;
    for (int ngen=0;ngen<=G;ngen++){
        // Nfg: first monomer's index of generation ngen
        // Ngl-1: last  monomer's index of generation ngen
        // ngen=0: Ngf=0;Ngl-1=1;
        // ngen=1: Ngf=2;Ngl-1=5;
        // ngen=2: Ngf=6;Ngl-1=13;
        Ngl   = Ngf + 2*(int) pow((F-1),(ngen));
        for (int m=Ngf;m<Ngl;m++){
            for (nbond=0;nbond<(F-1);nbond++){
                if (ngen!=G){
                    if (ngen == (G-1)) Mon[mnext].type = S;
                    else               Mon[mnext].type = C;
                    Mon[mnext].id        = mnext;
                    Mon[mnext].gen       = ngen+1;
                    Mon[mnext].bond[F-1] = m;
                    Mon[m].bond[nbond]   = mnext++;
                }
                else if (ngen == G)
                    Mon[m].bond[nbond] = -2;
            }
        }
        Ngf = Ngl;
    }
}

void InitConfig(){
    int    ok=0;
    int    countmon,m,i,mnext,nbond,ngen,Ngf,Ngl;
    double radius,tol,r2;
    int counter;
    VecR   vecr,dr;
    countmon=0;
    nbond = 0;
    ngen=0;
//
//  place gen=0 monomers, central monomers
    radius=fnBB.L0;
    // set the position of the first monomer to
    // the centre of the Box
    Mon[0].pos.x     =0.0;// (double) (Box.x/2);
    Mon[0].pos.y     =0.0;// (double) (Box.y/2);
    Mon[0].pos.z     =0.0; //(double) (Box.z/2);
    Mon[0].bond[F-1] = 1;
    Mon[0].id        = 0;
    Mon[0].type      = B;
    Mon[0].gen =0;
    RandomVector(&vecr,radius);
    Mon[1].pos.x     = Mon[0].pos.x + vecr.x;
    Mon[1].pos.y     = Mon[0].pos.y + vecr.y;
    Mon[1].pos.z     = Mon[0].pos.z + vecr.z;
    Mon[1].bond[F-1] = 0;
    Mon[1].id        = 1;
    Mon[1].type      = B;
    Mon[1].gen       = 0;
    

    Ngf   = 0;
    mnext = 2;  // monomers 0 and 1 are already placed
    for (ngen=0;ngen<G;ngen++){
        // Nfg: first monomer's index of generation ngen
        // Ngl-1: last  monomer's index of generation ngen
        // ngen=0: Ngf=0;Ngl-1=1;
        // ngen=1: Ngf=2;Ngl-1=1;
        Ngl   = Ngf + 2*(int) pow((F-1),(ngen));
        // loop over all monomers of generation ngen and
        // build monomers of generation ngen+1
        for (m=Ngf;m<Ngl;m++){
            if ((Mon[m].type==B || Mon[m].type==C)&&ngen==(G-1))
                                                radius = fnCS.L0;
            else                                radius = fnCC.L0;
            for (nbond=0;nbond<(F-1);nbond++){
                if (ngen == (G-1)) Mon[mnext].type = S;
                else               Mon[mnext].type = C;
                Mon[mnext].id        = mnext;
                Mon[mnext].gen       = ngen+1;
                Mon[m].bond[nbond]   = mnext;
                Mon[mnext].bond[F-1] = m;
                counter = 0;
                for (;;){
                    counter++;
                    RandomVector(&vecr,radius);
                    Mon[mnext].pos.x = Mon[m].pos.x+vecr.x;
                    Mon[mnext].pos.y = Mon[m].pos.y+vecr.y;
                    Mon[mnext].pos.z = Mon[m].pos.z+vecr.z;
                    ok = 1;
                    for (i=0;i<mnext;i++){
//                        printf("i=%-3d\n",i);
                        dr.x=Mon[mnext].pos.x-Mon[i].pos.x;
                        dr.y=Mon[mnext].pos.y-Mon[i].pos.y;
                        dr.z=Mon[mnext].pos.z-Mon[i].pos.z;
                        r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                        if (radius>1.0) tol = 1.0;
                        else            tol = SQR(radius);
                        if   (r2<tol)  {
                            ok = 0;
                            break;
                        }
                        else ok = 1;
                    }
                    if (ok == 1 ) break;
                }

                mnext++;
            }
        }
        Ngf  = Ngl;
    }
}

void RandomVector (VecR* vecr,double radius){
    // Marsaglia Algorithm
    // Random vector on a sphere
    //  of a certain radius
    // *** check  precision e.g r=1.000000045
    double ran1,ran2,ransq,ranh;
    
    for (;;){
        ran1=1.0-2.0*dr250();
        ran2=1.0-2.0*dr250();
        ransq=SQR(ran1)+SQR(ran2);
        if (ransq<1.0) break;
    }
    vecr->z=radius*(1.0-2.0*ransq);

    ranh=2.0*sqrt(1-ransq);
    radius=radius*ranh;
    
    vecr->x=radius*ran1;
    vecr->y=radius*ran2;

}
void TestIniConfig(FILE * pOut){
    for (int m=0;m<Nmon;m++){
        fprintf(pOut,"%d\t%d\t%8.3f\t%8.3f\t%8.3f\t\n",
                Mon[m].type,Mon[m].gen,Mon[m].pos.x,Mon[m].pos.y,Mon[m].pos.z);
    }
    fprintf(pOut,"\n Bonds atom: bond1, bond2, ....\n");
    for (int m=0;m<Nmon;m++){
        fprintf(pOut,"%3d: ",m);
        for (int b=0;b<F;b++){
            fprintf(pOut,"%3d "  "",Mon[m].bond[b]);
        }
        fprintf(pOut,"\n");
    }
}
void WriteXYZ (FILE *pFile){

    int dummy;

    dummy = ftell(pFile);
    if (dummy == 0)
        fprintf(pFile,"%-5d \n \n",Nmon);
    else
        fprintf(pFile,"\n\n");

    for (int m=0;m<Nmon;m++){
        fprintf(pFile,"G%-4d\t%8.3lf\t%8.3lf\t%8.3lf\n",
                Mon[m].gen,Mon[m].pos.x,Mon[m].pos.y,Mon[m].pos.z);
/*
        fprintf(pFile,"ATOM  %5d  %s %s  %4d  %8.3lf%8.3lf%8.3lf\n",
               m,"C  ",((Mon[m].type==0)?"B":((Mon[m].type==1)?"C":((Mon[m].type==2)?"S":"X"))),1,
                Mon[m].pos.x,Mon[m].pos.y,Mon[m].pos.z);
*/
    }
}
void WriteBondVmd (FILE* pFile){
//    fprintf(pFile,"mol new MC.xyz\n");
    fprintf(pFile,"mol new MC_G%dF%d.xyz\n",G,F);
    fprintf(pFile,"mol delrep 0 top\n");
    fprintf(pFile,"mol representation Bonds 0.1 10; mol addrep top\n");
    fprintf(pFile,"mol representation VDW   0.4 20; mol addrep top\n");
    fprintf(pFile,"set sel [ atomselect top \"all\" ]\n");
    fprintf(pFile,"$sel setbonds {\n");
    for (int mon=0;mon<Nmon;mon++){
        fprintf(pFile,"{ %3d ",mon);
        for (int nbond=0;nbond<F;nbond++){
            if (Mon[mon].bond[nbond]>=0){
                fprintf(pFile," %3d ",Mon[mon].bond[nbond]);
            }
        }
        fprintf(pFile,"}\n");
    }
    fprintf(pFile,"}\n");

}

int McMove (int NumofDispl){
    // move function for one Monte Carlo Step
    int  m,mtrial,accepted;
    VecR *oldPos,newPos;
    double energyOld,energyNew,dE;

    accepted =0;
    for(m=0;m<NumofDispl;m++){
        // choose random monomer
        mtrial = r250n(Nmon);
        // set old position as the position
        // of the selected monomer
        oldPos=&Mon[mtrial].pos;
        // calculate the energy of the monomer for
        // the old position
        EnergyMonomer(*oldPos,mtrial,0,&energyOld);
/*
        printf("Energy old: %25.3lf",energyOld);
*/
        TrialMove(&newPos,oldPos);

        EnergyMonomer(newPos,mtrial,0,&energyNew);
/*
        printf("Energy new: %25.3lf\n",energyNew);
*/
        dE=energyNew-energyOld;
        
        if (dE<0) {
            oldPos->x=newPos.x;
            oldPos->y=newPos.y;
            oldPos->z=newPos.z;
            RunningEnergy+=dE;
            accepted++;
        }

        else if (dr250()<exp(-dE)){
            oldPos->x=newPos.x;
            oldPos->y=newPos.y;
            oldPos->z=newPos.z;
            RunningEnergy+=dE;
            accepted++;
        }
/*
        else printf("Move Rejected: %80.3lf > %80.3lf \n",energyNew,energyOld);
*/
    }
    return (accepted);
}

void TrialMove (VecR *newPos, const VecR* oldPos){
    double radius;
    VecR   vecr;
    radius = MaxStep*dr250();

    RandomVector(&vecr,radius);

    newPos->x=oldPos->x+vecr.x;
    newPos->y=oldPos->y+vecr.y;
    newPos->z=oldPos->z+vecr.z;


}
void EnergyMonomer (VecR pos,int mcur,int mb,double*Ener){
// Calculates the energy of the mcur Monomer
// with monomers with id from mb to Nmon
// pos    : trial position of current monomer mcur
// typecur: type of current monomer Mon[mcur].type
    int    j,mnei,m,nbond,typecur,typenei;
    int    bonded=1;
    int    interType;
    VecR   dr;
    double en;
    double r,r2;
    FeneParams *pfn;
    MorseParams *pmrs;
    j=0;
    en=0;
    typecur = Mon[mcur].type;
    // FENE for bonded monomers
    for (nbond=0;nbond<F;nbond++){
        mnei    = Mon[mcur].bond[nbond];
        typenei = Mon[mnei].type;
        if (typecur == typenei){
            if      (typecur == B) pfn = &fnBB;
            else if (typecur == C) pfn = &fnCC;
        }
        if (typecur !=typenei){
            if (typecur == S || typenei == S) pfn = &fnCS;
            else                              pfn = &fnCC;
        }
        if ((mnei>=mb)&&(mnei!=-2)){
            dr.x=pos.x-Mon[mnei].pos.x;
            dr.y=pos.y-Mon[mnei].pos.y;
            dr.z=pos.z-Mon[mnei].pos.z;

            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            r=sqrt(r2);
            
            en+=Fene(pfn,r);
/*
            printf("mon1=%3d,mon2=%d\n",Mon[mcur].type,Mon[mnei].type);
            printf("potential param type->%s\n",pfn->type);
            printf("mon dist r=%lf\n",r);
            printf("params\n K=%lf, R=%lf, L0=%lf \n",pfn->K,pfn->R,pfn->L0);
            printf("energy en=%lf\n",en);
*/
            
/*
            en+=Fene(mcur,r);//-K*SQR(R)*log(1-SQR((r-L0)/R));
*/

        }
    }
    // Calculate non bonded monomer interactions that are
    // not already calculated
    for (m=mb;m<Nmon;m++){
        if (mcur!=m){
            for (nbond=0;nbond<F;nbond++){
                mnei    = Mon[mcur].bond[nbond];
                bonded = ((m==mnei) ? 1 : 0);
                if (bonded) break;
            }
            
            if (!bonded){
                typenei = Mon[m].type;
                if (typecur == typenei){
                    if (typecur == B || typecur == C) {
                        pmrs = &MrsCC;
#if USELUT==1
                        interType = CC;
#endif
                    }
                    if (typecur == S){
                        pmrs = &MrsSS;
#if USELUT==1
                        interType = SS;
#endif
                    }
                }
                if (typecur != typenei){
                    if (typecur == S || typenei == S) {
                        pmrs = &MrsCS;
#if USELUT==1
                        interType = CS;
#endif
                    }
                    else{
                        pmrs = &MrsCC;
#if USELUT==1
                        interType = CC;
#endif
                    }
                }
                dr.x=pos.x-Mon[m].pos.x;
                dr.y=pos.y-Mon[m].pos.y;
                dr.z=pos.z-Mon[m].pos.z;
                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
#if USELUT==1
                if (r2<TableRCutOff[interType]){
                double weight,rbin;
                int bin;
                rbin = (r2-TableRsqMin)*TableInvDeltaRsq[interType];
                bin  = (int)floor(rbin);
                if (bin <0) {
                    printf("negative bin number\n");
                    bin=0;
                }
                weight = rbin - bin;
                en += weight*MorseTab[interType][bin+1]+
                        (1-weight)*MorseTab[interType][bin];
                }
#else
                if (r2<SQR(RCut)){
                    r=sqrt(r2);
                    en +    = Morse(pmrs,r)-Morse(pmrs,RCut);
                }
#endif
            }
        }
    }
    *Ener=en;
}

double Fene(FeneParams *fnp,double r){
    return -(fnp->K) * SQR(fnp->R)*log(1-SQR((r-(fnp->L0))/(fnp->R)));
}

double Morse (MorseParams* mrsp,double r){
    return (mrsp->eps)*(SQR(exp(-(mrsp->a)*(r-(mrsp->d)))-1)-1);
}

void EnergySystem (){
    double en;
    en=0.0;
    
    EnergyTotal =0.0;

    for (int m=0;m<Nmon-1;m++){
        EnergyMonomer(Mon[m].pos,m,m+1,&en);
        EnergyTotal+=en;
    }
    // possible tail correction to the energy
    // should go here
}

void CalcProps (FILE *pProps,FILE *pCMass){
    static int propCounter=1;
    CalcCMass();
    WriteCMass(pCMass);
    CalcRg();
    fprintf(pProps,"%-10d\t%8.3f\t%8.3f\n",propCounter,RunningEnergy/Nmon,Rg.val);
    CalcDpNew();
    ++propCounter;


#if DEBUG ==1
    int m;
    FILE *pTest;
    pTest = fopen("cMass_Rg.dat","w");
    fprintf(pTest,"Positions of Monomers\n\n");
    for (m=0;m<Nmon;m++){
        fprintf(pTest,"pos  : %7.3lf  %7.3lf  %7.3lf\n",Mon[m].pos.x,Mon[m].pos.y,Mon[m].pos.z);
    }
    fprintf(pTest,"\n CMass: %7.3lf  %7.3lf  %7.3lf\n",CMass.x,CMass.y,CMass.z);
    fprintf(pTest,"\n Radius of gyration^2:%7.3lf \n",Rg2.val);
    fprintf(pTest,"\n Radius of gyration:%7.3lf \n",Rg.val);
#endif

}

void CalcCMass (){
    // Centre of Mass Calculation
    CMass.x=0.0;
    CMass.y=0.0;
    CMass.z=0.0;

    for (int m=0;m<Nmon;m++){
        CMass.x+=Mon[m].pos.x;
        CMass.y+=Mon[m].pos.y;
        CMass.z+=Mon[m].pos.z;
    }
    CMass.x/=(double)Nmon;
    CMass.y/=(double)Nmon;
    CMass.z/=(double)Nmon;
}

void WriteCMass (FILE *pCMass){
    static int step=1;
    fprintf(pCMass,"%-10d\n",step++);
    fprintf(pCMass,"%10.4f\t%10.4f\t%10.4f\n",CMass.x,CMass.y,CMass.z);
}
void CalcDpNew(){
    static int dpCount=0;
    VecR   dr;
    double rr,r,rg;
    int    m,h,ngen;
    double vb,normFac;
    if (dpCount == 0){
        for (ngen=0;ngen<=G;ngen++)
            for (h=0;h<dpSizeHist;h++)
                dpHistGen[ngen][h]=0.0;
    }
    deltaR=dpRange/dpSizeHist;
    for (m=0;m<Nmon;m++){
        ngen = Mon[m].gen;
        dr.x = Mon[m].pos.x-CMass.x;
        dr.y = Mon[m].pos.y-CMass.y;
        dr.z = Mon[m].pos.z-CMass.z;
        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if (rr<SQR(dpRange)){
            h=(int)floor(sqrt(rr)/deltaR);
                ++dpHistGen[ngen][h];
        }
    }
    ++dpCount;
    if (dpCount == dpLimit){
        // Calculations according to the average value of
        // radius of gyration
        rg = Rg.sum/Rg.count;
        for (ngen=0;ngen<=G;ngen++){
            for (h=0;h<dpSizeHist;h++){
                r=((double)h+0.5)*deltaR;
                vb = 4*PI*SQR(r)*deltaR;
                normFac=1.0/(vb*dpCount);
                dpHistGen[ngen][h]*=normFac;
            }
            WriteDpNew(ngen);
        }
        TestDpNew(stdout);
        dpCount = 0;
    }
        
    
}
void   WriteDpNew(int ngen){
    double r,rg;
    char filename[30];
    FILE  *pDp;
    if (ngen<=G){
        sprintf(filename,"dpHistG%dF%d.gen%d.dat",G,F,ngen);
        pDp = fopen(filename,"w");
        rg=Rg.sum/Rg.count;
        for (int h=0;h<dpSizeHist;h++){
            r=((double)h+0.5)*deltaR;
            fprintf(pDp,"%8.4f\t%8.4f\n",r/rg,dpHistGen[ngen][h]*CUBE(rg));
        }
        fclose(pDp);
    }
}
void TestDpNew (FILE* pFile){
    printf("radius of gyration:%8.3f\n",Rg.sum/((double)Rg.count));
    double r,vb;
    double dummy=0.0;
    for (int ngen=0;ngen<=G;ngen++)
        for (int h=0;h<dpSizeHist;h++){
            r=((double)h+0.5)*deltaR;
            vb = 4*PI*SQR(r)*deltaR;
            dummy+=dpHistGen[ngen][h]*vb;
        }
    fprintf(pFile,"**Testing Density Profile**\n");
    fprintf(pFile,"Total Num of Monomers DP:%-8.2f\n",dummy);
}
void CalcDp (){
    static int dpCount=0;
    VecR   dr;
    double rr,r,rg;
    int    m,h,t;
    double vb,normFac;

    if (dpCount == 0){
        for (t=0;t<NTYPE;t++)
            for (h=0;h<dpSizeHist;h++)
                dpHist[t][h]=0.0;
    }
    deltaR=dpRange/dpSizeHist;
    for (m=0;m<Nmon;m++){
        t    = Mon[m].type;
        dr.x = Mon[m].pos.x-CMass.x;
        dr.y = Mon[m].pos.y-CMass.y;
        dr.z = Mon[m].pos.z-CMass.z;

        rr=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
        if (rr<SQR(dpRange)){
            h=(int)floor(sqrt(rr)/deltaR);
            if (t == B || t == C ) {
                t=0;
                ++dpHist[t][h];
            }
            if (t == S){
                t=1;
                ++dpHist[t][h];
            }
/*
            if (m == 0)
    printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %3d\n",CMass.x,CMass.y,CMass.z,Mon[m].pos.x,Mon[m].pos.y,Mon[m].pos.z,h);
*/
        }
    }
    ++ dpCount;

    // Normalization
    if (dpCount == dpLimit){
        // Calculations according to the average value of
        // radius of gyration
        rg = Rg.sum/Rg.count;
        for (t=0;t<NTYPE;t++){
            for (h=0;h<dpSizeHist;h++){
                r=((double)h+0.5)*deltaR;
                vb = 4*PI*SQR(r)*deltaR;
                normFac=1.0/(vb*dpCount);
                dpHist[t][h]*=normFac;
            }
        }
        WriteDp(0);
        WriteDp(1);
        TestDp(stdout);
        dpCount = 0;
    }
}

void TestDp (FILE* pFile){
    printf("radius of gyration:%8.3f\n",Rg.sum/((double)Rg.count));
    double r,vb;
    double dummy=0.0;
    for (int t=0;t<NTYPE;t++)
        for (int h=0;h<dpSizeHist;h++){
            r=((double)h+0.5)*deltaR;
            vb = 4*PI*SQR(r)*deltaR;
            dummy+=dpHist[t][h]*vb;
        }
    fprintf(pFile,"**Testing Density Profile**\n");
    fprintf(pFile,"Total Num of Monomers DP:%-8.2f\n",dummy);
}

void WriteDp(int t){
    double r,rg;
    char filename[30];
    FILE  *pDp;
    if (t<NTYPE){
        sprintf(filename,"dpHist%c.dat",(t==0 ? 'C' : (t==1 ? 'S':'N')));
        pDp = fopen(filename,"w");
        rg=Rg.sum/Rg.count;
        for (int h=0;h<dpSizeHist;h++){
            r=((double)h+0.5)*deltaR;
            fprintf(pDp,"%8.4f\t%8.4f\n",r/rg,dpHist[t][h]*CUBE(rg));
        }
        fclose(pDp);
    }
}

void CalcRg (){
    static int dpRg=0;
    int m;

    if (dpRg == 0) {
        Rg.sum    = 0.0;
        Rg.count  = 0;
        Rg2.sum   = 0.0;
        Rg2.count = 0;
    }

    Rg2.val = 0.0;
    for (m=0;m<Nmon;m++){
        Rg2.val+=SQR(Mon[m].pos.x-CMass.x)+
                 SQR(Mon[m].pos.y-CMass.y)+
                 SQR(Mon[m].pos.z-CMass.z);
    }
    ++dpRg;
    Rg2.val /= (double) Nmon;
    Rg2.sum += Rg2.val;
    ++Rg2.count;
/*
    printf("Rg:%8.3f\t%-3d",Rg.sum,Rg.count);
*/
    Rg.val   = sqrt(Rg2.val);
    Rg.sum  += Rg.val;
    ++Rg.count;
}

void AdjustStep(double accRatio){
    static int att=0,acc=0;
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
        if (MaxStep/curstep>1.5) MaxStep = curstep*1.5;
        if (MaxStep/curstep<0.5) MaxStep = curstep*0.5;
//        if (MaxStep/curstep>Box.x/4.0) MaxStep = Box.x/4.0;
#if DEBUG == 1
        printf("Old Max Step: %lf\n",curstep);
        printf("New Max Step: %lf\n",MaxStep);
        printf("Adjusting acc. ratio to:%lf\n",accRatio);
        printf("Attempts:%d\n",NumofAttempts-att);
        printf("Accepted:%d\n",NumofAcceptedMoves-acc);
#endif

    }

}

void BuildCellList (void){

    int n,c,cyz,numofCells;
    cells.x=Box.x/RCut;
    cells.y=Box.y/RCut;
    cells.y=Box.y/RCut;

    cyz=cells.y*cells.z;
    numofCells=cells.x*cells.y*cells.z;

    rc.x=Box.x/cells.x;
    rc.y=Box.y/cells.y;
    rc.z=Box.z/cells.z;

    cellList = (int*) malloc((numofCells*Nmon)*sizeof(int));

    for (c=Nmon;c<Nmon+numofCells;c++) cellList[c]=EMPTY;
    for (n=0;n<Nmon;n++) {
        cInd.x=Mon[n].pos.x/rc.x;
        cInd.y=Mon[n].pos.y/rc.y;
        cInd.z=Mon[n].pos.z/rc.z;

        c=cInd.x*cyz+cInd.x*cells.z+cInd.z+Nmon;

        cellList[n]=cellList[c];
        cellList[c]=n;
    }

}


