 /* 
 * File:   system.h
 * Author: georgiou
 *
 * Created on February 9, 2011, 4:07 PM
 */

#ifndef SYSTEM_H
#define	SYSTEM_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define NTYPE  2   // number of different monomer types (core, shell)
#define B      0   // central monomers
#define C      1   // core monomers 
#define S      2   // shell monomers 

#define PI 3.141592653589793

typedef struct{int x;} VecTest;
typedef struct {int x, y, z;}  VecI;
typedef struct {double x,y,z;} VecR;
typedef struct {
    double val,sum;
    unsigned long int count;
} dbProp;
// define true false for bool variables
#define VCopy(v1, v2)                                       \
   (v1).x = (v2).x,                                         \
   (v1).y = (v2).y,                                         \
   (v1).z = (v2).z
#define VSub(v1, v2, v3)                                    \
   (v1).x = (v2).x - (v3).x,                                \
   (v1).y = (v2).y - (v3).y,                                \
   (v1).z = (v2).z - (v3).z
#define VDot(v1, v2)                                        \
   ((v1).x * (v2).x + (v1).y * (v2).y + (v1).z * (v2).z)
#define VLenSq(v)  VDot (v, v)
#define VLen(v)  sqrt (VDot (v, v))
#define TRUE (1==1)
#define FALSE (!TRUE)

#define SQR(x)  ( (x) * (x) )
#define CUBE(x) ( (x) * (x) * (x) )
#define MIN(x1, x2)  (((x1) < (x2)) ? (x1) : (x2))
#define MAX(x1, x2)  (((x1) > (x2)) ? (x1) : (x2))
#define MAX3(x1,x2,x3)   \
   (((x1) > (x2)) ? (((x1) > (x3)) ? (x1) : (x3)) :         \
                    (((x2) > (x3)) ? (x2) : (x3)))
#define AllocMat(mat,rows,cols,t) mat = (t *) malloc(rows * cols * sizeof(t))
#define A(mat,rows,cols,i,j) mat[i*cols + j]
	

// Periodic Boundary Conditions
// RAPAPORT definitions
#define VWrap(v, t)                                   \
   if (v.t >= 0.5 * Box.t)      v.t -= Box.t;         \
   else if (v.t < -0.5 * Box.t) v.t += Box.t

#define VWrapAll(v)                                   \
   {VWrap (v, x);                                     \
   VWrap (v, y);                                      \
   VWrap (v, z);}

#define PBC1(v,t) ( (v).t = (v).t - rint(((v).t)/Box.t) * (Box.t) )

#define PBCAll(v)   \
{   PBC1(v,x);      \
    PBC1(v,y);      \
    PBC1(v,z); }
#define PBC_(v0,v,t) ( (v0).t = (v).t - rint(((v).t)/Box.t)*(Box.t) )
// Variables in system.c
//
extern char   fileIniConf[120];
extern char   filePotential[120];
extern int    Restart;
extern int    TaskID;
extern unsigned int  seed;
extern FILE  **pOutputFile; // one file for each dendrimer
extern unsigned int  numOfFramesFile;

extern int    G,F,numOfXYZfiles; //generation(G),spacer length(P) not needed here,
								 //functionality(F)
extern double iniDendDist;
extern VecI   Box;
extern double CellListSizeMaxrad;
extern double Density;
extern double Temp, TempHigh,TempLow,TempStep;
extern double Beta;

extern double EnergyTotal;
extern double MorseTotal, FeneTotal;
extern double RunningEnergy;
extern double MaxStep;
extern int    PotType;

extern double *RsqCut; // Holds RCut Values for Morse Potential
                     //Rcut[N_TYPE_MORSE] :: CC->0 CS->1 SS->2

// dendrimers.c
//*****************************************//
//  Dendrimer - Monomers definittions      //
struct strMonomer{
    int    id;
    int    uniqid;
    int    type;
    int    inDend;
    int    gen;
    int    numOfBonds;
    double rad; // radius for vmd
	VecR   trialpos;
    VecR   pos;
    VecR   npos; // pointer uncorrected positions
                 // if no PBC are used it points to pos
    VecR   bupos;
    // Bonds & Bond Interaction
    struct strMonomer **bond;
    int    *bondInterType;
    // linked cell
    int cellIndex; // position in hoc
    struct strMonomer *next;
    struct strMonomer *previous;
};

typedef struct strMonomer Monomer;

struct strDendrimer{
    int     id;
    int     gen;
    int     func;
    //int     first; // index of first monomer in dendrimer
    //int     last;  // indes of last !!! for (mon=first;mon <= last; mon++)
    VecR    CMass;
    dbProp  Rg;
    dbProp  Rg2;
	double  b,c,k2,Rgsq;
    int     numOfMonomers;
    int     *numOfMonGen;
    struct strMonomer *monomer;
    double  rNext;  // distance of next monomer
};
typedef struct strDendrimer Dendrimer;

extern int numOfDendrimers;
extern int numOfMonomersInDend;
extern int numOfTotalMonomers;
extern int numOfTotalBonds;
extern int numOfIniConfigFiles;

//extern  Monomer   *monomer;
extern  Dendrimer *dendrimer;

// dendrimers.c
//*****************************************//
//             Dendrimers                  //
//*****************************************//
void    DendSetup(void);
int 	DendReadCoords(const char inFileName[],int);
int     DendSaveConfig(const char outFileName[],int);
void    DendDefineBondsTypes(Dendrimer* );
void    DendCalcCMass(Dendrimer* );
void    DendCalcRg(Dendrimer* );
double *DendGyrTensor(Dendrimer*);
void    DendInitConfig(int , const VecR*);
void    DendInitConfig2(int , const VecR*);
void    DendWriteProps(const Dendrimer*, FILE * );
double  DendMaxRad(const Dendrimer *pd );
void    DendWriteXYZ(void);
void    DendWriteXYZCommon(FILE*,int);
void    DendWriteCMassXYZ(FILE*,int);
void    WriteBondVmd (const char[]);
void  	DendWriteinlmp(const char [],const char [],const int );
void    DendWriteBonds (FILE*);
int     DendCheckBox(void);
void    DendUnfoldXYZ(void);
void    DendPBCunwrap(VecR*, const VecR*);
void    DendSetnpos(Dendrimer*);
void    DendSetpos(Dendrimer*);
int     IsEqual(double, double, double);
void    EnergyCheck(void);

//potential.c definitions
//*****************************************//
//           Potential Parameter           //
//*****************************************//
typedef struct{
    int    interType;
    char   type[3];
    double K;
    double L0;
    double R;
}FeneParams;


typedef struct{
    int    interType;
    char   type[3];
    double eps;
    double a;
    double d;
    double RCut;
    double RLow;
}MorseParams;

#define N_TYPE_MORSE     3      // number of different interaction types
#define N_TYPE_FENE      3      // for MORSE and FENE Potentials
//**********Interaction type**************//
// use the same symbols in verlet lists
// core - core interaction for both
// FENE and MORSE Potential
#define CC 0
#define CS 1
#define SS 2 //
#define BB 2 // !!!!!**use only in FENE**!!!!!

extern FeneParams  fnCC,fnCS,fnBB;
extern MorseParams MrsCC,MrsCS,MrsSS;

double Fene(FeneParams* ,double);     //exact calculation of FENE potential
double Morse (MorseParams* ,double ); //exact calculation of Morse potential
void   PotentialSetup(int);            // the seven Mladek Potentials

void   MorseOutputValues(FILE*);
void   WriteFene(int, double, double, double, const char[]);
void   readtoarray(const char [],double [],double [],int*);
double readfromarray(double ,double [],double [],const int);
//-----------------------------------------//
//energy.c definitions
//*****************************************//
double InterDend(Dendrimer*, Dendrimer*);
double InterDendRef(Dendrimer, Dendrimer);
double InterDendId(int, int);
double IntraDend(Dendrimer*);
double EnergySystem(void);// Calculates total energy of system of
                          //interacting dendrimers
double MorseEnergy(Monomer*, VecR*);
double FeneEnergy(Monomer*,VecR*,int);
double PotMonIntraDend(Monomer*, VecR * ,int);
double PotMonNextDend(Monomer*, VecR * , Dendrimer*);
double PotDendDend(Dendrimer *, Dendrimer*);
double EnergySystem2(void);
double EnergySystem3(void);
double EnergyDend (Dendrimer*);
double EnergyMonomer (VecR pos,Monomer *pm,int mb);
void CalcInterForce(Dendrimer* , Dendrimer* ,VecR *,VecR *);
double  MorseForce(MorseParams* ,VecR );

double CelldESingleMonMove(Monomer* , VecR *);
void ApplyBoundaryConditions(void);

//-----------------------------------------//
//moves.c definitions
//*****************************************//
void DendPlaceCMass(Dendrimer*, VecR);
void DendMoveCMass(Dendrimer *,VecR);
void MonTrialMoveFixedCMass(Monomer *);
void MonSingleMove (Monomer *, VecR*, int);
void MonTrialMove(VecR*, const VecR*); // option to keep centre of mass fixed
void MonTrialMove_All(VecR *, const VecR *,VecR* , const VecR *);
void RandomVector(VecR*,double);
void WriteVecR(VecR, FILE*);
void WriteVecI(VecI, FILE*);


//-----------------------------------------//
//*****************************************//
//                 mc.c                    //
#define DECORRELATION  0
#define EQUILIBRATION  1
#define PRODUCTION     2

extern long int NumofDecorrelationCycles;
extern long int NumofEquilibrationCycles;
extern long int NumofProductionCycles;
extern int      FramesFrequency;
extern long int NumofAttempts;
extern long int NumofAcceptedMoves;
extern long int NumofCycles;
extern int      NumofDispPerCycle;
extern double   AcceptanceRatio;
extern int      SamplingFrequency;

void MC(void);
void MCanneal(double ,double , double );
void MCequil(int , int );
void MC_multi(void);
int  McMove(int);
int  McMove_two(int);
int  McMove_t(int);
void AdjustStep(const double);
void MCPropEnergy(FILE*);
void MCProps(FILE*, FILE*, FILE*, FILE*);

//**************************************//
//         Functions definitions        //
//**************************************//
void ReadSimParams(const char []);
void ReadSimParamsAnneal(const char []);
void ReadPotential(int ,const char []);
void OutputInputParams(void);
void Initialize(void);
void SetParams(void);
void Allocate(void);

//**************************************//
//             Cell Lists               //
//**************************************//
extern double  maxRc; // maximum cutoff of all Morse potential
struct strCell{
    int index;
    int numOfNeigh;
    struct strMonomer *first;
    int numOfElements;
    int neighbours[125];
};

typedef struct strCell Cell;
extern  struct strCell *cells;
extern VecI    cellNI;                   // number of cells in directions x,y,z;
extern int     cellNYZ;                  // cellNx * cellNy
extern VecR    cellSizeR, cellInvSizeR;  // size and inverse size of cell in each direction
extern int     cellTotalNumber;          // Total number of cells in the box

#define IFLOAT(a,L) ( (int) floor( (a)/(L)) )
#define DCELLINDX(ax,ay,az,CS) CELLINDX(IFLOAT(ax,(CS).x),IFLOAT(ay,(CS).y),IFLOAT(az,(CS).z),cellNI)
// index = cx + cy * Lx + cz*Lx*Ly
#define CELLINDX(xx,yy,zz,CN) ( (xx) + (yy)*((CN).x) + (zz)*((CN).x)*((CN).y) )
// cx = mod (index,Lx), cy=mod[(index/ Lx),Ly] ,cz = index/(Lx*Ly)
#define CELLCOORDS(VI,cindx,CN)                       \
{   VI.x =    (cindx) % (CN).x;                       \
    VI.y =  ( (cindx) / ((CN).x) ) % ((CN).y);        \
    VI.z =    (cindx) / ( ((CN).y) * ((CN).z) );}     \

#define VCELLINDX(a) DCELLINDX((a).x,(a).y,(a).z,CS)

void CellListSetup(VecI,double);
void CellAlloc(void);
void CellDeAlloc(void);
void CellNewBuild(void);
void CellUpdateMon(struct strMonomer*);
void CellUpdateAll(struct strDendrimer*);
int  CellUpdateList();// update all dendrimers
void CellAddElement(int, struct strMonomer*);
void CellDelElement(int, struct strMonomer*);
void CellTest(void);
double CellInterDend(void);
void CellOutElements(int);
// LUT DEFINITIONS
#define LUT_MORSE_SIZE 10000
#define LUT_FENE_SIZE  10000
extern double **MorseLUT,*MLUTRsqMin,*MLUTRsqCut,*MLUTDeltaRsq,
			   *MLUTInvDeltaRsq,*MLUTRCut,*MLUTPhiCutoff;

void   LUTMorseAllocate(void); // allocate look-up tables tables
void   LUTMorseBuild(int);     // calculate values of LUT
void   LUTMorseOutputSetup(const char[]);
void   LUTMorseOutputTables(const char[],double,double,double);
double LUTMorseValue(int ,double); // input: interaction type ,squared distance,
                                           // pointer to double energy
// Fene Look-up
extern double **FeneLUT,*FLUTRsqMin,*FLUTRsqMax,*FLUTDeltaRsq,*FLUTInvDeltaRsq,
               *FLUTPhiRsqMin,*FLUTPhiRsqMax;

void   LUTFeneAllocate(void);
void   LUTFeneBuild(double);     // calculate values of LUT
void   LUTFeneOutputTables(const char[],double,double,double);
void   LUTFeneOutputSetup(const char[]);
double LUTFeneValue(int ,double);
#endif	/* SYSTEM_H */
