#include "system.h"

char fileIniConf[120];
char filePotential[120];
int  Restart;
int  TaskID;
unsigned  int seed;

VecI   Box;
double CellListSizeMaxrad;
double Density;
double Temp, TempHigh,TempLow,TempStep;
double Beta;
double *RsqCut;
    
int    G,F;
double iniDendDist;
double EnergyTotal;
double MorseTotal, FeneTotal;
double RunningEnergy;
double MaxStep;
int    PotType;

int    numOfDendrimers;
int    numOfMonomersInDend;
int    numOfTotalMonomers;
int	   numOfTotalBonds;
int    numOfXYZfiles;
int    numOfIniConfigFiles;

//struct strMonomer *monomer;
struct strDendrimer *dendrimer;
// Cell list Neighbouhr
//*******************************//
//              MC               //
long int NumofDecorrelationCycles;
long int NumofEquilibrationCycles;
long int NumofProductionCycles;
int      FramesFrequency;
long int NumofAttempts;
long int NumofAcceptedMoves;
long int NumofCycles;
int      NumofDispPerCycle;
double   AcceptanceRatio;
int      SamplingFrequency;

FeneParams  fnCC,fnCS,fnBB;
MorseParams MrsCC,MrsCS,MrsSS;

double **MorseLUT = NULL,*MLUTRCut = NULL,*MLUTRsqMin = NULL,
        *MLUTDeltaRsq = NULL,*MLUTInvDeltaRsq = NULL,*MLUTRsqCut = NULL,
        *MLUTPhiCutoff = NULL;

double **FeneLUT = NULL, *FLUTRsqMin = NULL, *FLUTRsqMax = NULL ,
        *FLUTDeltaRsq = NULL ,*FLUTInvDeltaRsq = NULL,*FLUTPhiRsqMin = NULL,
        *FLUTPhiRsqMax = NULL;

//cellLists,c
struct strCell *cells=NULL;
//Monomer **hoc=NULL;

double  maxRc;
VecI    cellNI;
int     cellNYZ;
VecR    cellSizeR, cellInvSizeR;
int     cellTotalNumber;
