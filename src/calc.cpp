/* 
 * Author: georgiou
 * Created on August 11, 2011, 1:36 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>
#include <cmath>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "nr3.h"
#include "eigen_sym.h"
#include "ran250.h"

#define CUBE(x)  ((x)*(x)*(x))
#define MSQR(x) ((x)*(x))
#define _PI 3.1415926535897932384626433832
#define TEST 0
#define FABS(x) (x)
#include "rotations.h"

typedef struct{
    double x,y,z;
}VecR;

typedef struct{
    int    id;
    VecR   CMass;
    double Rg,b,c,k2;
    double RgC,bC,cC,k2C;
    double RgS,bS,cS,k2S;
	int P2;
    VecR   laxis;
	vector<VecR> moncoords;
	vector<VecR> monforce;// total force acting to monomer
	double forceCM; // force on CMass (for 2 interacting dendrimers only)
	VecR   vforceCM;
	VecR   vforceCMtot;
}Dend;

void PotentialSetup(int);
double dot_product(const VecR& va,const VecR& vb){
    double result =   (va.x)*(vb.x)+(va.y)*(vb.y)+(va.z)*(vb.z);
    return result;
};

using namespace std;
int    numOfDendrimers, G, F, pottype;

long long int frames;

vector< vector<int> >  numofMonGen;
vector<int>numofCoreMons (numOfDendrimers,0);
vector<int>numofShellMons(numOfDendrimers,0);

vector<double> sumb,sumbC,sumbS ;
vector<double> sumsqb,sumsqbC,sumsqbS;
vector<double> sumc,sumcC,sumcS;
vector<double> sumsqc,sumsqcC,sumsqcS;
vector<double> sumRg,sumsqRg,averRg;
vector<double> sumRgC,sumsqRgC,averRgC;
vector<double> sumRgS,sumsqRgS,averRgS;

VecR   CMass;
double   Rg;

void CalcEigen(const vector<VecR>& pvcoords,int mysize,VecDoub&d,MatDoub&v);
void CalcCMass(const vector<VecR> &pvcoords,VecR & rCMass);
void CalcRg(const vector<VecR>&,const VecR* );
void CalcMsdCM(VecR& ,VecR& ,double &);
void CalcMsdMon(const vector<VecR>& ,const VecR& ,const vector<VecR>& ,const VecR& ,double &);
void CalcDp(const int& ,const double, const VecR& ,const vector<VecR> &, vector< vector< vector<double> > > &);
void CalcP2(vector<double>& pvP2,vector< vector<int> >& order,vector<Dend>& myD,int numOfPairs);
void ScaleDP(const int& , const int &,const double ,vector< vector< vector<double> > > &);
void TestDp (double , vector< vector< vector<double> > > &);

// FORCE CALCULATION DECLARONS

#define B      0   // central monomers
#define C      1   // core monomers 1st generation
#define S      2   // shell monomers second generation

#define CC 0
#define CS 1
#define SS 2 //
#define BB 2 // !!!!!**use only in FENE**!!!!!

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

FeneParams  fnCC,fnCS,fnBB;
MorseParams MrsCC,MrsCS,MrsSS;
void MorseForce(MorseParams& pmrs,VecR& dr);
void FeneForce(VecR& ff,VecR& dr);
double CalcForce(Dend& d1, Dend& d2);
int  Montype(vector<int> & pvcormon,vector<int>& pvshellmon,int dendid,int monid);
double testTotalForce(vector<Dend> & D);
// END OF FORCE CALCULATION DECLARATIONS

void outputCoords(const vector<VecR>& ,const char *);
void RandomSpherePoint(vector<VecR>& ,int,double );
void RandomRotateCoords(vector<VecR>& );
void RandomVector(VecR& , double ) ;

bool OpenCppFileExists(const string& filename);


int main(int argc, char** argv) {
    if (argc<9){
        cout<<"too few arguments"<<endl;
        cout<<"(program_name) [numofdendrimers] [G] [F] [pottype]"
                "[ini_file_name in xyz format)] [output b c file] [distance] [num of frames] [dpHistSize] [dpRange])"<<endl;
        exit(0);
    }
//  initilize random number generator 
    long int seed = 30982;
    srand(seed);
    // Input file contains particles positions x, y, z
    // the program reads only the first three columns
    std::string   infilename, outfilename;
    std::ifstream infile;
    
    // Output Files
    std::ofstream outrun,outbc, outCMass, outProps, outEigen, outOrient, outMSD, outP2,outForceFrame;
    stringstream filenameCMass, filenameProps, filenameEigen, filenameOrient, filenameMSD,filenameP2,filenameForceFrame;
    std::string   sline;

    int    numOfMonomersInDend,numOfTotalMonomers;
    double dummyX, dummyY, dummyZ;
    double dpRange;
    
//  Gyration tensor stuff
    const  int gyrtensorsize=3;
//  Shape of whole molecule
    VecDoub d(gyrtensorsize); // 1D array returning the eigen values of matrix
    MatDoub v(gyrtensorsize,gyrtensorsize);// 2D array where columns are the eigenvectors
//  Shape of Core and Shell seperately    
    VecDoub dC(gyrtensorsize),dS(gyrtensorsize);
	MatDoub vC(gyrtensorsize,gyrtensorsize),vS(gyrtensorsize,gyrtensorsize);
    
    double  iniDist; // only for interactions between two dendrimers
    VecR    centralbond, largestAxis;// central bond vector ,coordinates of largest axis
    
    double Rgsq, normadot, b, c, k2, angle;
    double RgsqC, bC, cC, k2C,angleC;
    double RgsqS, bS, cS, k2S,angleS;
    // dot products
    double dot_central_laxis,dot_central_fixed,dot_laxis_fixed;
    int dpHistSize = 1000;

    // Orientation parameters
    // Cmass for dendrimers 0 and 1
    // Vectors connecting the two CMasses starting from 0 or 1 respectively
    VecR CMass0, CMass1, vCMass0, vCMass1;
    // angle between vCMass and long axis of dendrimers 0 or 1
    // angle between long axes of dendrimers 0 and 1


    // counters
    int dendid, mon, gen, countmonomers;
	long long int countframes;
    // Read arguments store to variablesvoid CalcMsdCM(VecR& rcurCMass,VecR& rrefCMass,double &rdblMsdCM){
    
    int ii=1;
    
    numOfDendrimers = atoi(argv[ii]);
    G               = atoi(argv[++ii]);
    F               = atoi(argv[++ii]);
	pottype         = atoi(argv[++ii]);
    infilename      = string(argv[++ii]);
    outfilename     = string(argv[++ii]);
    iniDist         = atof(argv[++ii]);
    frames          = atoi(argv[++ii]);
    dpHistSize      = atoi(argv[++ii]);
    dpRange         = atof(argv[++ii]);
    outrun.open("run.out",std::ios::app);
	outrun<<"processing file:"<<infilename<<endl;
	outrun<<"NumDend:"<<numOfDendrimers<<endl
          <<"TypeD"<<pottype<<endl         
          <<"frame:"<<frames<<endl         
          <<"histsize"<<dpHistSize<<endl
          <<"histrange"<<dpRange<<endl<<endl;         
	outrun.close();

    numOfMonomersInDend = 2*(1-(int)pow((F-1),G+1))/(2-F);
#if TEST==1
    numOfMonomersInDend = 6;
#endif
    numOfTotalMonomers  = numOfDendrimers * numOfMonomersInDend;
	PotentialSetup(pottype);

    // alocate vectors
    // numofMonGen
    vector<Dend>  D(numOfDendrimers);
	for (int dd=0;dd<D.size();dd++){
		D.at(dd).monforce.resize(numOfMonomersInDend);
	}	
    vector<VecR> moncoords(numOfMonomersInDend);
    vector< vector<VecR> > refcoords(numOfDendrimers,vector<VecR> (numOfMonomersInDend));

    // Density Profile
    // [dendrimer id][generation][histogram size]
    // create 3D vector of double to hold histogram
    vector< vector< vector<double> > > dpHistDendGen(numOfDendrimers, vector< vector<double> >((G+1), vector<double>(dpHistSize, 0)));
    // mean square displacement vector for each dendrimer
    vector<double> msd(numOfDendrimers,0.0);
    //vector< vector<double> > msdmon(numOfDendrimers,vector<double>(numOfMonomersInDend,0.0)); 
    vector<double> msdCMass(numOfDendrimers,0.0);
    vector<double> msdMon(numOfDendrimers,0.0);
    const int numOfDendPairs=numOfDendrimers*(numOfDendrimers-1)/2;;
    vector<double> P2(numOfDendPairs,0.0);    
	vector< vector<int> > orderPair(numOfDendPairs,vector<int>(2,0));// order[pair][0]=dend0 of pair ; order[pair][1]=dend1 of pair; initialized in CalcP2
    // CMass reference vector
    vector<VecR> refCMass(numOfDendrimers);

	vector< vector<VecR> > dendmoncoords( numOfDendrimers, vector<VecR> (numOfMonomersInDend) );
    
    // specify number of monomers in each generation

    numofMonGen.resize(numOfDendrimers);

	numofCoreMons.resize(numOfDendrimers); 
	numofShellMons.resize(numOfDendrimers); 
    for (dendid =0;dendid < numOfDendrimers;dendid++){
        numofMonGen.at(dendid).resize(G+1);
		D.at(dendid).vforceCM.x=0.0;
		D.at(dendid).vforceCM.y=0.0;
		D.at(dendid).vforceCM.z=0.0;
		D.at(dendid).vforceCMtot.x=0.0;
		D.at(dendid).vforceCMtot.y=0.0;
		D.at(dendid).vforceCMtot.z=0.0;
        for (gen = 0; gen <= G; gen++) {
                numofMonGen.at(dendid).at(gen) =2 * (int) pow((F-1), gen);
                if (gen<G)  numofCoreMons.at(dendid) += 2 * (int) pow((F-1), gen);
                if (gen==G) numofShellMons.at(dendid)+= 2 * (int) pow((F-1), gen);
        }
    }

	// projected force from all frames
	VecR F01={0.0,0.0,0.0};
	VecR F10={0.0,0.0,0.0};
	// force on CM for current frame
	VecR F01fr={0.0,0.0,0.0};
	VecR F10fr={0.0,0.0,0.0};
    // asphericity and acylindricity
      sumb.resize(numOfDendrimers,0.0);
    sumsqb.resize(numOfDendrimers,0.0);
      sumc.resize(numOfDendrimers,0.0);
    sumsqc.resize(numOfDendrimers,0.0);
    
      sumbC.resize(numOfDendrimers,0.0);
    sumsqbC.resize(numOfDendrimers,0.0);
      sumcC.resize(numOfDendrimers,0.0);
    sumsqcC.resize(numOfDendrimers,0.0);
    
      sumbS.resize(numOfDendrimers,0.0);
    sumsqbS.resize(numOfDendrimers,0.0);
      sumcS.resize(numOfDendrimers,0.0);
    sumsqcS.resize(numOfDendrimers,0.0);

	  sumRg.resize(numOfDendrimers,0.0);
	sumsqRg.resize(numOfDendrimers,0.0);
	 averRg.resize(numOfDendrimers,0.0);

	  sumRgC.resize(numOfDendrimers,0.0);
	sumsqRgC.resize(numOfDendrimers,0.0);
	 averRgC.resize(numOfDendrimers,0.0);

	  sumRgS.resize(numOfDendrimers,0.0);
	sumsqRgS.resize(numOfDendrimers,0.0);
	 averRgS.resize(numOfDendrimers,0.0);
    
    infile.open(infilename.c_str(),std::ios::in);
    if (infile.fail()){
        cout<<"cannot open input file"<<infilename.c_str()<<endl;
        exit(0);
    }
    
    filenameCMass<<"CMassn"<<numOfDendrimers<<"G"<<G<<"F"<<F<<"d"<<iniDist<<".dat";
    filenameEigen<<"Eigenn"<<numOfDendrimers<<"G"<<G<<"F"<<F<<"d"<<iniDist<<".dat";
    filenameProps<<"Propsn"<<numOfDendrimers<<"G"<<G<<"F"<<F<<"d"<<iniDist<<".dat";
    filenameOrient<<"Orientn"<<numOfDendrimers<<"G"<<G<<"F"<<F<<"d"<<iniDist<<".dat";
    filenameMSD<<"MSDn"<<numOfDendrimers<<"G"<<G<<"F"<<F<<"d"<<iniDist<<".dat";
    filenameP2<<"P2cos"<<numOfDendrimers<<"G"<<G<<"F"<<F<<"d"<<iniDist<<".dat";
	filenameForceFrame<<"F_fr"<<"G"<<G<<"F"<<F<<"d"<<iniDist<<".dat";
    
    outbc.open(outfilename.c_str(),std::ios::app);
    
    outCMass.open(filenameCMass.str().c_str(),std::ios::out);
    outEigen.open(filenameEigen.str().c_str(),std::ios::out);
    outProps.open(filenameProps.str().c_str(),std::ios::out);
    outOrient.open(filenameOrient.str().c_str(),std::ios::out);
    outMSD.open(filenameMSD.str().c_str(),std::ios::out);
    outP2.open(filenameP2.str().c_str(),std::ios::out);
	outForceFrame.open(filenameForceFrame.str().c_str(),std::ios::out);

    outCMass<<"#"<<endl<<"#"<<endl<<"#"<<endl;

    outProps<<noshowpos<<"# num_dend"<<"/"<<"num_mon"<<"/"
            <<"G"<<"/"<<"F"<<"/"<<"num of frames read"<<endl;
    outProps<<noshowpos<<"#"<<numOfDendrimers<<"/"<<numOfMonomersInDend<<"/"
            <<G<<"/"<<F<<"/"<<frames<<endl;
    outProps<<"#  1:id   2:b/Rgsq    3:c/Rgsq       4:k2      5:Rgsq       6:Rg    7:angle(laxis,central)"; // All monomers
    outProps<<"#  8:b/Rgsq(C)  9:c/Rgsq (C)   10:k2(C)  11:Rgsq(C)    12:Rg(C) ";// Core  monomers
    outProps<<"# 14:b/Rgsq(S)  15:c/Rgsq(S)   16:k2(S)  17:Rgsq(S)    18:Rg(S)    ";// Shell monomers
    
            
    outEigen<<noshowpos<<"# num_dend"<<"/"<<"num_mon"<<"/"
            <<"G"<<"/"<<"F"<<"/"<<"num of frames read"<<endl;
    outEigen<<noshowpos<<"#"<<numOfDendrimers<<"/"<<numOfMonomersInDend<<"/"<<
            G<<"/"<<F<<"/"<<frames<<endl;
    outEigen<<"#D_id\t"<<"    X^2           Y^2          Z^2\t"<<
              "                             eX                "<<
              "                             eY                "<<
              "                             eZ                "<<endl;
	outP2<<"di dj P2[cos(theta)]"<<endl;
	outForceFrame<<"1:frame 2:F01_x 3:F01_y 4:F01_z 5:F10_x 6:F10_y 7:F10_z "<<endl;
    countframes = 0;
    int readfreq=1; // read ever readfreqframe
#if TEST==1
    // create random points on a sphere
        RandomSpherePoint(moncoords,numOfMonomersInDend,1.0);
        while(countframes<frames){
#else
        while ((infile.good())&&(countframes<frames)&&(!infile.eof())) {
#endif
        for (dendid = 0;dendid<numOfDendrimers;dendid++){
            countmonomers=0;
#if TEST==1
            //rotate coords using quaternions
            RandomRotateCoords(moncoords);
            countmonomers+=numOfMonomersInDend;
            outputCoords(moncoords,"test.out");
#else            
                for (mon=0;mon<moncoords.size();mon++){
                    getline(infile,sline);
                    // using old good c
                    sscanf(sline.c_str(), "%lf %lf %lf", &dummyX, &dummyY, &dummyZ);
                    moncoords.at(mon).x = dummyX;
                    moncoords.at(mon).y = dummyY;
                    moncoords.at(mon).z = dummyZ;                
                    countmonomers++;
                }
			//	dendmoncoords.at(dendid)=moncoords;
				D.at(dendid).moncoords=moncoords;				
				moncoords=D.at(dendid).moncoords;
#endif            
			//	moncoords=D.at(dendid).moncoords;
             // store reference coords for Min square displacement calculations
                if (countframes==0){
                
                for (mon=0;mon<refcoords.at(dendid).size();mon++){
                    refcoords.at(dendid).at(mon).x = moncoords.at(mon).x;
                    refcoords.at(dendid).at(mon).y = moncoords.at(mon).y;
                    refcoords.at(dendid).at(mon).z = moncoords.at(mon).z;
                }
                CalcCMass(refcoords.at(dendid),refCMass.at(dendid));
            }
            if ((countframes%readfreq)==0){

//# START CALCULATIONS
            CalcCMass(moncoords,CMass);
            CalcRg(moncoords,&CMass);
            CalcMsdCM(CMass,refCMass.at(dendid),msdCMass.at(dendid));
            CalcMsdMon(moncoords,CMass,refcoords.at(dendid),refCMass.at(dendid),msdMon.at(dendid));
            
            CalcEigen(moncoords,numOfMonomersInDend,d,v);
            CalcEigen(moncoords,numofCoreMons.at(dendid),dC,vC);
            CalcEigen(moncoords,numofShellMons.at(dendid),dS,vS);
            
            Rgsq = d[0] + d[1] + d[2];
            b    = d[0] - 0.5*( d[1] + d[2] );
            c    = d[1] - d[2];
            k2   = ( MSQR(b) + 0.75 * MSQR(c) )/ MSQR(Rgsq);
            
              sumb.at(dendid) += b/Rgsq;
            sumsqb.at(dendid) += MSQR(b/Rgsq);
              sumc.at(dendid) += c/Rgsq;
            sumsqc.at(dendid) += MSQR(c/Rgsq);

			  sumRg.at(dendid) += sqrt(Rgsq);
			sumsqRg.at(dendid) += (Rgsq);
            
            RgsqC = dC[0] + dC[1] + dC[2];
            bC    = dC[0] - 0.5*( dC[1] + dC[2] );
            cC    = dC[1] - dC[2];
            k2C   = ( MSQR(bC) + 0.75 * MSQR(cC) )/ MSQR(RgsqC);
            
              sumbC.at(dendid) += bC/RgsqC;
            sumsqbC.at(dendid) += MSQR(bC/RgsqC);
              sumcC.at(dendid) += cC/RgsqC;
            sumsqcC.at(dendid) += MSQR(cC/RgsqC);

			  sumRgC.at(dendid) += sqrt(RgsqC);
			sumsqRgC.at(dendid) += (RgsqC);
            
            RgsqS = dS[0] + dS[1] + dS[2];
            bS    = dS[0] - 0.5*( dS[1] + dS[2] );
            cS    = dS[1] - dS[2];
            k2S   = ( MSQR(bS) + 0.75 * MSQR(cS) )/ MSQR(RgsqS);
            
              sumbS.at(dendid) += bS/Rgsq;
            sumsqbS.at(dendid) += MSQR(bS/Rgsq);
              sumcS.at(dendid) += cS/Rgsq;
            sumsqcS.at(dendid) += MSQR(cS/Rgsq);

			  sumRgS.at(dendid) += sqrt(RgsqS);
			sumsqRgS.at(dendid) += (RgsqS);
            
            D.at(dendid).CMass.x = CMass.x;
            D.at(dendid).CMass.y = CMass.y;
            D.at(dendid).CMass.z = CMass.z;
            D.at(dendid).Rg      = sqrt(Rgsq);
            D.at(dendid).b       = b;
            D.at(dendid).c       = c;
            D.at(dendid).k2      = k2;
            D.at(dendid).RgC     = sqrt(RgsqC);
            D.at(dendid).bC      = bC;
            D.at(dendid).cC      = cC;
            D.at(dendid).k2C     = k2C;
            D.at(dendid).RgS     = sqrt(RgsqS);
            D.at(dendid).bS      = bS;
            D.at(dendid).cS      = cS;
            D.at(dendid).k2S     = k2S;
            D.at(dendid).laxis.x = v[0][0];
            D.at(dendid).laxis.y = v[1][0];
            D.at(dendid).laxis.z = v[2][0];
            
//# END CALCULATIONS            
            if (numOfDendrimers==2){
                if (dendid==0){
                    CMass0.x = CMass.x;
                    CMass0.y = CMass.y;
                    CMass0.z = CMass.z;

                }
                if (dendid==1){
                    CMass1.x = CMass.x;
                    CMass1.y = CMass.y;
                    CMass1.z = CMass.z;
                    // calculate the vector connecting the
                    // two centre of Masses
                    vCMass0.x = CMass1.x - CMass0.x;
                    vCMass0.y = CMass1.y - CMass0.y;
                    vCMass0.z = CMass1.z - CMass0.z;
                    
                    vCMass1.x = CMass0.x - CMass1.x;
                    vCMass1.y = CMass0.y - CMass1.y;
                    vCMass1.z = CMass0.z - CMass1.z;

                    double normadotvCMass0, normadotvCMass1;
                    normadotvCMass0 = sqrt(FABS(dot_product(vCMass0,vCMass0)));
                    normadotvCMass1 = sqrt(FABS(dot_product(vCMass1,vCMass1)));
                    
                    vCMass0.x /= normadotvCMass0;
                    vCMass0.y /= normadotvCMass0;
                    vCMass0.z /= normadotvCMass0;

                    vCMass1.x /= normadotvCMass1;
                    vCMass1.y /= normadotvCMass1;
                    vCMass1.z /= normadotvCMass1;
                    
                }
            }
            outCMass<<setw(4)<<setprecision(0)<<noshowpoint<<noshowpos<<fixed
                    <<countframes+1<<" "<<dendid<<"\t";
            outCMass<<showpoint<<showpos<<scientific<<setw(13)<<setprecision(6)
                    <<CMass.x<<"  "<<CMass.y<<"  "<<CMass.z<<" "<<sqrt(MSQR(CMass.x)+MSQR(CMass.y)+MSQR(CMass.z))<<endl;
            
            
            
            outProps<<setw(4)<<setprecision(0)<<noshowpoint<<noshowpos<<fixed<<dendid<<" ";
            outProps<<showpoint<<scientific<<setw(11)<<setprecision(5)
                    <<b/Rgsq<<" "<<c/Rgsq<<" "<<k2<<" "<<Rgsq<<" "<<sqrt(Rgsq)<<" "<<angle<<endl;
                    //<<bC/RgsqC<<" "<<cC/RgsqC<<" "<<k2C<<" "<<RgsqC<<" "<<sqrt(RgsqC)<<" "
                    //<<bS/RgsqS<<" "<<cS/RgsqS<<" "<<k2S<<" "<<RgsqS<<" "<<sqrt(RgsqS)<<endl;
// END OF EIGENVALUES CALCULATIONS
            CalcDp(dendid,dpRange, D.at(dendid).CMass,D.at(dendid).moncoords, dpHistDendGen);
            
            centralbond.x = moncoords.at(1).x - moncoords.at(0).x;
            centralbond.y = moncoords.at(1).y - moncoords.at(0).y;
            centralbond.z = moncoords.at(1).z - moncoords.at(0).z;

            largestAxis.x = v[0][0];
            largestAxis.y = v[1][0];
            largestAxis.z = v[2][0];

//          Calculate angle between central bond and largest axis

            normadot = sqrt(FABS(dot_product(centralbond,centralbond)));
            centralbond.x /= normadot;
            centralbond.y /= normadot;
            centralbond.z /= normadot;
            angle = acos(FABS(dot_central_laxis)) * 180 / M_PI ;
            if (numOfDendrimers==1){
                VecR vfixed={1.0,0.0,0.0};
                if (dendid==0){
                    largestAxis.x = v[0][0];
                    largestAxis.y = v[1][0];
                    largestAxis.z = v[2][0];
                    // all vectors have to be normalized
                    // central bond with long axis
                    dot_central_laxis = dot_product(largestAxis,centralbond);
                    // central bond with const vector vfixed
                    dot_central_fixed = dot_product(centralbond,vfixed);
                    // long axis with const vector vfixed
                    dot_laxis_fixed   = dot_product(largestAxis,vfixed);
                    double angle_central_laxis = acos( FABS(dot_central_laxis) ) * 180 / M_PI ;
                    double angle_central_fixed = acos( FABS(dot_central_fixed) ) * 180 / M_PI ;
                    double angle_laxis_fixed   = acos( FABS(dot_laxis_fixed)   ) * 180 / M_PI ;

                    if ((countframes==0)&&(numOfDendrimers==1)){
                        outOrient<<"#frame  2:angle(central,laxis) "
                                   "3:angle(central,fixed) "
                                   "4:angle(laxis,fixed) "
                                   "5,6,7:centralbond(xyz) "
                                   "8,9,10:laxis(x,y,z)"
                                 <<endl;
                    }
                    
                    
                    outOrient.setf(ios::fixed);
                    outOrient<<setw(5)<<setprecision(0)<<noshowpoint<<noshowpos
                            <<countframes+1<<"\t";
                    outOrient<<showpoint<<setw(9)<<setprecision(2)
                         <<angle_central_laxis<<"\t"
                         <<angle_central_fixed<<"\t"
                         <<angle_laxis_fixed<<"\t"
                         <<setw(6)<<setprecision(4)
                         <<centralbond.x<<"\t"<<centralbond.y<<"\t"<<centralbond.z<<"\t"
                         <<largestAxis.x<<"\t"<<largestAxis.y<<"\t"<<largestAxis.z<<endl;
                }
            }
            if (numOfDendrimers==2){
                VecR ex0, ex1;            
                VecR centralB0,centralB1;
                if (dendid==0){
                    ex0.x = v[0][0];
                    ex0.y = v[1][0];
                    ex0.z = v[2][0];
                    centralB0.x = moncoords.at(1).x - moncoords.at(0).x;
                    centralB0.y = moncoords.at(1).y - moncoords.at(0).y;
                    centralB0.z = moncoords.at(1).z - moncoords.at(0).z;
                    normadot = sqrt(FABS(dot_product(centralB0,centralB0)));
                    centralB0.x /= normadot;
                    centralB0.y /= normadot;
                    centralB0.z /= normadot;
                }
                if (dendid==1){
                    ex1.x = v[0][0];
                    ex1.y = v[1][0];
                    ex1.z = v[2][0];
                    centralB1.x = moncoords.at(1).x - moncoords.at(0).x;
                    centralB1.y = moncoords.at(1).y - moncoords.at(0).y;
                    centralB1.z = moncoords.at(1).z - moncoords.at(0).z;
                    normadot = sqrt(FABS(dot_product(centralB1,centralB1)));
                    centralB1.x /= normadot;
                    centralB1.y /= normadot;
                    centralB1.z /= normadot;
                    // all vectors have to be normalized
                    double angle_laxis01       = acos( FABS(dot_product(ex0,ex1)) )     *180 /M_PI;
					double dummydot = MSQR(dot_product(ex0,ex1));
                    double angle_laxis0_CMass  = acos( FABS(dot_product(vCMass0,ex0)) ) *180 /M_PI;
                    double angle_laxis1_CMass  = acos( FABS(dot_product(vCMass1,ex1)) ) *180 /M_PI;
                    
                    VecR dr;double rsq;
                    dr.x = vCMass0.x;
                    dr.y = vCMass0.y;
                    dr.z = vCMass0.z;
                    rsq = MSQR(dr.x)+MSQR(dr.y)+MSQR(dr.z);
                    //output
                    if ((countframes==0)){
                        outOrient<<"#frame  2:angle(laxis0,laxis1) "
                                   "3:angle(laxis0,CMassCon01) "
                                   "4:angle(laxis1,CMassCon10) "
                                   "5,6,7    : centralbond0(xyz)"
                                   "8,9,10   : centralbond1(xyz)"
                                   "11,12,13 : laxis0(x,y,z)"
                                   "14,15,16 : laxis1(x,y,z)"
                                   "17: angledif[angle(laxis0,CMassCon01),angle(laxis1,CMassCon10)] "
								   "18: (laxis0*laxis1)^2"
                                
                                 <<endl;
                    }
                    outOrient.setf(ios::fixed);
                    outOrient<<right;
                    outOrient<<setw(5)<<setprecision(0)<<noshowpoint<<noshowpos<<countframes+1<<"\t";
                    outOrient<<showpoint<<showpos<<setw(8)<<setprecision(2)
                             <<angle_laxis01<<"\t"<<angle_laxis0_CMass<<"\t"<<angle_laxis1_CMass<<"\t"
                             <<centralB0.x<<"\t"<<centralB0.y<<"\t"<<centralB0.z<<"\t"
                             <<centralB1.x<<"\t"<<centralB1.y<<"\t"<<centralB1.z<<"\t"
                             <<ex0.x<<"\t"<<ex0.y<<"\t"<<ex0.z<<"\t"
                             <<ex1.x<<"\t"<<ex1.y<<"\t"<<ex1.z<<"\t"
                             <<fabs(angle_laxis0_CMass-angle_laxis1_CMass)<<dummydot<<endl;
                }
                
            }

            
            
            
            // OUTPUT FORMAT
            // For each dendrimer
            // d[0] d[1] d[2] v[i][0] v[i][1] v[i][2] i=0..2
            // end of frames to empty lines
            outEigen<<setw(3)<<setprecision(0)<<noshowpoint<<noshowpos<<fixed<<dendid<<"\t";
            
            for (int i=0;i<d.size();i++)
                outEigen<<showpoint<<showpos<<fixed<<setw(13)<<setprecision(6)<<d[i]<<" ";
            outEigen<<"\t";
            double norm;
            for (int col=0;col<v.ncols();col++){
                norm=0;
                for (int row=0;row<v.nrows();row++){
                    norm+=MSQR(v[row][col]);
                    outEigen<<showpoint<<showpos<<scientific<<setw(13)<<setprecision(6)<<v[row][col]<<" ";
                }
                //outEigen<<norm<<" "<<"\t";
                outEigen<<"\t";
            }
            outEigen<<endl;
                
        }
        }
        CalcP2(P2,orderPair,D,numOfDendPairs);
		if (numOfDendrimers==2)
		{


			VecR vCM01,vCM10;
			double norma01,norma10;
			
			vCM01.x = D.at(0).CMass.x - D.at(1).CMass.x;
			vCM01.y = D.at(0).CMass.y - D.at(1).CMass.y;
			vCM01.z = D.at(0).CMass.z - D.at(1).CMass.z;
			vCM01.x = vCM01.x/sqrt(dot_product(vCM01,vCM01));
			vCM01.y = vCM01.y/sqrt(dot_product(vCM01,vCM01));
			vCM01.z = vCM01.z/sqrt(dot_product(vCM01,vCM01)); 	

			double rr =MSQR(vCM01.x)+MSQR(vCM01.y)+MSQR(vCM01.z);
			if (iniDist< 0.4){
					vCM01.x =1.0;			
					vCM01.y =0.0;
					vCM01.z =0.0;
			}
			// unit vector

			vCM10.x = -vCM01.x;
			vCM10.y = -vCM01.y;
			vCM10.z = -vCM01.z;

//			printf("%lf %lf %lf\n",	vCM10.x,vCM10.y,vCM10.z);
			CalcForce(D.at(0),D.at(1));		

			for (int did=0;did<numOfDendrimers;did++){
					D.at(did).vforceCM.x = 0.0;
					D.at(did).vforceCM.y = 0.0;
					D.at(did).vforceCM.z = 0.0;
					for (mon=0;mon<D.at(0).moncoords.size();mon++){
							D.at(did).vforceCMtot.x += D.at(did).monforce.at(mon).x;
							D.at(did).vforceCMtot.y += D.at(did).monforce.at(mon).y;
							D.at(did).vforceCMtot.z += D.at(did).monforce.at(mon).z;
							D.at(did).vforceCM.x    += D.at(did).monforce.at(mon).x;
							D.at(did).vforceCM.y    += D.at(did).monforce.at(mon).y;
							D.at(did).vforceCM.z    += D.at(did).monforce.at(mon).z;
					}
			}

			double dummyProj;
			dummyProj = dot_product(D.at(0).vforceCM,vCM01)/dot_product(vCM01,vCM01);

			F01.x   += dummyProj*vCM01.x; 
			F01.y   += dummyProj*vCM01.y; 
			F01.z   += dummyProj*vCM01.z; 

			F01fr.x  = D.at(0).vforceCM.x; 
			F01fr.y  = D.at(0).vforceCM.y; 
			F01fr.z  = D.at(0).vforceCM.z; 
						
			dummyProj = dot_product(D.at(1).vforceCM,vCM10)/dot_product(vCM10,vCM10);
			F10.x += dummyProj*vCM10.x; 
			F10.y += dummyProj*vCM10.y; 
			F10.z += dummyProj*vCM10.z; 

			F10fr.x  = D.at(1).vforceCM.x; 
			F10fr.y  = D.at(1).vforceCM.y; 
			F10fr.z  = D.at(1).vforceCM.z; 
			/*D.at(0).forceCM +=sqrt(dot_product(F01,F01));
			D.at(1).forceCM +=sqrt(dot_product(F10,F10));
			D.at(0).vforceCM.x += F00.x;
			D.at(0).vforceCM.y += F01.y;
			D.at(0).vforceCM.z += F01.z;
			D.at(1).vforceCM.x += F10.x;
			D.at(1).vforceCM.y += F10.y;
			D.at(1).vforceCM.z += F10.z;*/
			//outForce<<iniDist<<" "<<sqrt(dot_product(F01,F01))<<" "<<sqrt(dot_product(F10,F10))<<" "<<F01.x<<" "<<F01.y<<" "<<F01.z<<" "<<F10.x<<" "<<F10.y<<" "<<F10.z<<endl;
			//outForce.close();
		}	
		//output force for current frame
		if (numOfDendrimers==2){
			outForceFrame<<setw(8)<<countframes+1<<" ";
			outForceFrame<<showpoint<<fixed<<showpos<<setw(15)<<setprecision(6); // force precission
			outForceFrame<<F01fr.x<<" "<<F01fr.y<<" "<<F01fr.z<<" "<<F10fr.x<<" "<<F10fr.y<<" "<<F10fr.z<<" "<<endl;
		}
        if (countframes==0) outMSD<<"#frame     dendid    MSD(dendid)"<<endl;
        for (dendid=0;dendid<numOfDendrimers;dendid++){
            outMSD<<setw(10)<<countframes+1<<" "<<setw(3)<<dendid<<" "
                  <<setw(10)<<setprecision(5)<<sqrt(msdCMass.at(dendid))<<" "
                  <<setw(10)<<setprecision(5)<<sqrt(msdMon.at(dendid))<<endl;
        }
        
        if (numOfDendrimers>1){
            outEigen<<endl<<endl;
            outProps<<endl<<endl;
            outCMass<<endl<<endl;
        }
        ++countframes;
        infile.ignore(256,'\n');
        infile.ignore(256,'\n');
        

        if (countmonomers!=numOfMonomersInDend){
            cout<<"Total monomers read("<<countmonomers<<")"<<"!= Total monomers in dendrimer("<<numOfMonomersInDend<<")"<<endl;
        }
    } // END OF FILE WHILE READ LOOP
	for (dendid=0;dendid<numOfDendrimers;dendid++){
		averRg.at(dendid) = sumRg.at(dendid)/countframes;	
	}
    outbc.seekp(0,ios::end);            
    if (outbc.tellp()==0){
        outbc<<"#1:dist  2:dendid   3:b  4:devb  5:c  6:devc 7:Rg 8:devRg "
               " 9:b(C) 10:devb(C) 11:c(C)  12:devc(C) 13:Rg(C)  14:devRg(C) "
               "15:b(S) 16:devb(S) 17:c(S)  18:devc(S) 19:Rg(S)  20:devRg(S) "<<endl;
    }        
	if (numOfDendrimers==2){
		//Project 
		stringstream filenameP2dist;
		ofstream outP2dist;
		filenameP2dist<<"P2vsDist.out";
		stringstream filenameForce;
		ofstream outForce;
		filenameForce<<"Force.out";
		outForce.open(filenameForce.str().c_str(),std::ios::app);
    	outForce.seekp(0,ios::end);            
    	if (outForce.tellp()==0){
    	    outForce<<"#1:dist  2:dist/Rg 3:F01 4:F10 "<<endl;
    	}        
		outP2dist.open(filenameP2dist.str().c_str(),std::ios::app);
		outP2dist<<iniDist<<" "<<iniDist/((D.at(0).Rg+D.at(1).Rg)/2.0)<<" "<<P2.at(0)/countmonomers<<" "<<0.5*(3*P2.at(0)/countframes-1.0)<<endl;
		/*outForce<<showpoint<<fixed<<showpos<<setw(13)<<setprecision(6)
				<<F01.x/countframes<<" "<<F01.y/countframes<<" "<<F01.z/countframes<<" "<<F10.x/countframes<<" "<<F10.y/countframes<<" "<<F10.z/countframes<<" "
				<<" "<<D.at(0).vforceCMtot.x/countframes<<" "<<D.at(0).vforceCMtot.y/countframes<<" "<<D.at(0).vforceCMtot.z/countframes
				<<" "<<D.at(1).vforceCMtot.x/countframes<<" "<<D.at(1).vforceCMtot.y/countframes<<" "<<D.at(1).vforceCMtot.z/countframes
				<<endl;*/
		outForce<<showpoint<<showpos<<setw(6)<<setprecision(3)<<iniDist<<" "<<iniDist/((D.at(0).Rg+D.at(1).Rg)/2.0)<<" ";
		outForce<<showpoint<<fixed<<showpos<<setw(13)<<setprecision(6)
				<<F01.x/countframes<<" "<<F10.x/countframes<<" "
				<<F01.y/countframes<<" "<<F10.y/countframes<<" "
				<<F01.z/countframes<<" "<<F10.z/countframes<<" "
				<<averRg.at(0)<<" "<<averRg.at(1)<<" "
				<<endl;
	}
   	double tempP2;
	for (int p2index=0;p2index<P2.size();p2index++){
		tempP2 = 0.5*(3*P2.at(p2index)/countframes-1);
		outP2<<p2index<<" "<<orderPair.at(p2index).at(0)<<" "<<orderPair.at(p2index).at(1)<<" "<<tempP2<<" "<<P2.at(p2index)<<endl;
	}
    double devb,devc,devbC,devcC,devbS,devcS,devRg,devRgC,devRgS;
    for (dendid=0;dendid<numOfDendrimers;dendid++){

		
        ScaleDP(dendid,countframes,dpRange,dpHistDendGen);
        if (iniDist>=0.0){
            devb  = sqrt( (sumsqb.at(dendid)/countframes) - MSQR(sumb.at(dendid)/countframes) );
            devc  = sqrt( (sumsqc.at(dendid)/countframes) - MSQR(sumc.at(dendid)/countframes) );
			devRg = sqrt( (sumsqRg.at(dendid)/countframes) - MSQR(sumRg.at(dendid)/countframes) );
            
            devbC = sqrt( (sumsqbC.at(dendid)/countframes) - MSQR(sumbC.at(dendid)/countframes) );
            devcC = sqrt( (sumsqcC.at(dendid)/countframes) - MSQR(sumcC.at(dendid)/countframes) );
			devRgC = sqrt( (sumsqRgC.at(dendid)/countframes) - MSQR(sumRgC.at(dendid)/countframes) );
            
            devbS = sqrt( (sumsqbS.at(dendid)/countframes) - MSQR(sumbS.at(dendid)/countframes) );
            devcS = sqrt( (sumsqcS.at(dendid)/countframes) - MSQR(sumcS.at(dendid)/countframes) );
			devRgS = sqrt( (sumsqRgS.at(dendid)/countframes) - MSQR(sumRgS.at(dendid)/countframes) );

            outbc<<iniDist<<" "<<dendid<<" "
                 <<sumb.at(dendid)/countframes<<" "<<devb<<" "<<sumc.at(dendid)/countframes<<" "<<devc<<" "<<sumRg.at(dendid)/countframes<<" "<<devRg<<" "
                 <<sumbC.at(dendid)/countframes<<" "<<devbC<<" "<<sumcC.at(dendid)/countframes<<" "<<devcC<<" "<<sumRgC.at(dendid)/countframes<<" "<<devRgC<<" "
                 <<sumbS.at(dendid)/countframes<<" "<<devbS<<" "<<sumcS.at(dendid)/countframes<<" "<<devcS<<" "<<sumRgS.at(dendid)/countframes<<" "<<devRgS<<" "<<endl;
        }
    }
			
			
    cout<<"Total Number of frames read:"<<countframes<<endl;

    outEigen.close();
    outCMass.close();
    outProps.close();
    outbc.close();
    outMSD.close();
	outP2.close();
    return 0;
}

void RandomSpherePoint(vector<VecR>& rvcoords,int numofmon,double radius){
    if (numofmon==6){
        //x direction 
            rvcoords.at(0).x =  1.0; rvcoords.at(0).y =  0.0;rvcoords.at(0).z =  0.0;
            rvcoords.at(1).x = -1.0; rvcoords.at(1).y =  0.0;rvcoords.at(1).z =  0.0;
            rvcoords.at(2).x =  0.0; rvcoords.at(2).y =  1.0;rvcoords.at(2).z =  0.0;
            rvcoords.at(3).x =  0.0; rvcoords.at(3).y = -1.0;rvcoords.at(3).z =  0.0;
            rvcoords.at(4).x =  0.0; rvcoords.at(4).y =  0.0;rvcoords.at(4).z =  1.0;   
            rvcoords.at(5).x =  0.0; rvcoords.at(5).y =  0.0;rvcoords.at(5).z = -1.0;            
    }else{
        for (int mon=0;mon<rvcoords.size();mon++){
            RandomVector(rvcoords.at(mon),radius);
        }
    }
}
void RandomRotateCoords(vector<VecR>& rvcoords){
    RMat   rotmat;
    Quat   quatRot;
    double eAng[3];
    // create a random eulerian angle
    RandomEulerAngle(eAng);
    //fprintf(stdout,"%lf %lf %lf\n",eAng[0],eAng[1],eAng[2]);
    // convert the angle to a quaternion
    EulerToQuat(&quatRot,eAng);
    // build rotation matrix
    BuildRotMatrix(&rotmat,&quatRot,1);
    // rotate vector positions of all monomers in molecule
    for (int mon=0;mon<rvcoords.size();mon++){
        VecR *pvecr,tmpvcr;
        pvecr = &rvcoords.at(mon);
        MVMul(tmpvcr,rotmat.u,*pvecr);
        rvcoords.at(mon).x = tmpvcr.x;
        rvcoords.at(mon).y = tmpvcr.y;
        rvcoords.at(mon).z = tmpvcr.z;
    }    
    
}
void RandomVector(VecR& rvecr, double radius) {
    // Marsaglia Algorithm
    // Random vector on a sphere
    //  of a certain radius
    // *** check  precision e.g r=1.000000045
    double ran1, ran2, ransq, ranh;

    for (;;) {
        ran1 = 1.0 - 2.0 * dr250();
        ran2 = 1.0 - 2.0 * dr250();
        ransq = SQR(ran1) + SQR(ran2);
        if (ransq < 1.0) break;
    }
    rvecr.z = radius * (1.0 - 2.0 * ransq);

    ranh = 2.0 * sqrt(1 - ransq);
    radius = radius*ranh;

    rvecr.x = radius*ran1;
    rvecr.y = radius*ran2;

}

void outputCoords(const vector<VecR>& rvcoords,const char *charfilename){
    static long int count=0;
    std::ofstream out;
    double radius;
    out.open(charfilename,std::ios::app);
    if (count>0){
        out<<endl<<endl;
    }
    for (int mon=0;mon<rvcoords.size();mon++){
        
        out<<setprecision(3)<<setw(8)<<showpos<<showpoint
            <<rvcoords.at(mon).x<<" "
            <<rvcoords.at(mon).y<<" "
            <<rvcoords.at(mon).z<<" ";
        radius = MSQR(rvcoords.at(mon).x)+MSQR(rvcoords.at(mon).y)+
                 MSQR(rvcoords.at(mon).z);
        //out<<radius<<endl;
        out<<0.5<<" "<<1<<" "<<1<<endl;;
    }
    count++;
}
void CalcDp(const int& dendid,const double dprange, const VecR& CMass,const vector<VecR> &rvcoords, vector< vector< vector<double> > > &hist){
    static int count=0;
    int  m, dpsize, mfirst, mlast, h;
    VecR dr;
    double rsq;
    dpsize = hist.at(0).at(0).size();
    double deltaR=dprange/dpsize;
    
    mfirst = 0; mlast = 0;
    for (int gen =0; gen<numofMonGen.at(dendid).size(); gen++){
        mlast = mfirst + numofMonGen.at(dendid).at(gen);
        for (m=mfirst;m<mlast;m++){
            dr.x = rvcoords.at(m).x-CMass.x;
            dr.y = rvcoords.at(m).y-CMass.y;
            dr.z = rvcoords.at(m).z-CMass.z;
            rsq=MSQR(dr.x)+MSQR(dr.y)+MSQR(dr.z);
            if (rsq<MSQR(dprange)){
                h=(int)floor(sqrt(rsq)/deltaR);
                hist.at(dendid).at(gen).at(h)++;
            }
            else{
                cout<<"!!!!!!!!! CalcDP()::distance out of range"<<endl;
            }
        }
        mfirst = mlast;
    }
    ++count;
}
void ScaleDP(const int& dendid, const int &frames,const double dprange,vector< vector< vector<double> > > &hist){
    double r,vb,normFac;
    int gen;
    double myconst, rupper,rlower; // variables for volume calculation
    stringstream histfilename;
    ofstream     outHist;
    int    Nbincor = 0; // number of bins to correct

    histfilename<<"hist_DP_d"<<dendid+1<<"g"<<G<<"f"<<F<<".dat";
    outHist.open(histfilename.str().c_str(),std::ios::out);

    int dpsize = hist.at(0).at(0).size();
    double deltaR = dprange/dpsize;

    myconst = 4.0*_PI/3.0;
    double dummyhist;
    // for the first [Nbincor] bins

    rupper = Nbincor*deltaR;
    rlower = 0;
    r =(rupper+rlower)/2.0;
    vb = myconst*(pow<double>(rupper,3.0)-pow<double>(rlower,3.0));
    
    normFac = 1.0/(vb*frames);

    outHist<<setw(13)<<setprecision(6)<<scientific<<r;
    
    for (gen=0;gen<numofMonGen.at(dendid).size();gen++){
        dummyhist=0.0;
        for (int h=0;h<Nbincor;h++){
            dummyhist+=hist.at(dendid).at(gen).at(h);
        }
        outHist<<" "<<dummyhist*normFac;
    }
    outHist<<endl;
    

    for (int h=Nbincor;h<dpsize;h++){
        rupper = (h+1)*deltaR;
        rlower = h*deltaR;
        r       = ((double)(h+0.5))*deltaR;
        vb      = 4*_PI*MSQR(r)*deltaR;
        //vb = myconst*(pow<double>(rupper,3.0)-pow<double>(rlower,3.0));
        normFac = 1.0/(vb*frames);
        outHist<<setw(13)<<setprecision(6)<<scientific<<r<<" "<<vb;
        for (gen=0;gen<numofMonGen.at(dendid).size();gen++){
//            hist.at(dendid).at(gen).at(h)*=normFac;
            outHist<<" "<<hist.at(dendid).at(gen).at(h)*normFac;
        }
        outHist<<endl;
    }

//    for (gen=0;gen<numofMonGen.at(dendid).size();gen++){
//        for (int h=0;h<dpsize;h++){
//            r=((double)h+0.5)*deltaR;
//            vb = 4*_PI*MSQR(r)*deltaR;
//            normFac=1.0/(vb*frames);
//            hist.at(dendid).at(gen).at(h)*=normFac;
//            outHist<<setw(13)<<setprecision(6)<<scientific<<r<<"\t"<<hist.at(dendid).at(gen).at(h)<<endl;
//        }
//        outHist<<endl<<endl;
//    }
    outHist.close();
}
void TestDp (double dprange, vector< vector< vector<double> > > &hist){
    double r,vb;
    double dummy=0.0;
    int    gen,dendid;
    int dpsize = hist.at(0).at(0).size();
    double deltaR = dprange/dpsize;
    for (int h=0;h<dpsize;h++){
        r=((double)h+0.5)*deltaR;
        vb = 4*_PI*MSQR(r)*deltaR;
        for (dendid=0;dendid<hist.size();dendid++){
            for (gen=0;gen<numofMonGen.at(dendid).size();gen++){
                dummy+=hist.at(dendid).at(gen).at(h)*vb;
            }
        }
    }
    cout<<"**Testing Density Profile**"<<endl;
    cout<<"Total Num of Monomers DP:%-8.2f\n"<<dummy<<endl;
}

void CalcCMass(const vector<VecR> &pvcoords,VecR& rCMass){
    int mon;
    // Centre of Mass Calculation
    rCMass.x = 0.0;
    rCMass.y = 0.0;
    rCMass.z = 0.0;

    for (mon = 0; mon < pvcoords.size(); mon++) {
        rCMass.x += pvcoords.at(mon).x;
        rCMass.y += pvcoords.at(mon).y;
        rCMass.z += pvcoords.at(mon).z;
    }
    rCMass.x /= (double) pvcoords.size();
    rCMass.y /= (double) pvcoords.size();
    rCMass.z /= (double) pvcoords.size();
}
void CalcRg(const vector<VecR>& pvcoords,const VecR* pCMass){
    int mon;
    Rg=0.0;
    double rg2=0.0;
    
    for (mon = 0; mon < pvcoords.size(); mon++) {
        rg2 += MSQR(pvcoords.at(mon).x-pCMass->x)+
               MSQR(pvcoords.at(mon).y-pCMass->y)+
               MSQR(pvcoords.at(mon).z-pCMass->z);
    }
    rg2 /= pvcoords.size();
    Rg= sqrt(rg2);
}
void CalcMsdCM(VecR& rcurCMass,VecR& rrefCMass,double &rdblMsdCM){
    VecR dr;
    dr.x = rcurCMass.x - rrefCMass.x;
    dr.y = rcurCMass.y - rrefCMass.y;
    dr.z = rcurCMass.z - rrefCMass.z;
    
    rdblMsdCM = MSQR(dr.x)+MSQR(dr.y)+MSQR(dr.z);
}
void CalcMsdMon(const vector<VecR>& pvcoords,const VecR& rCMass, 
                const vector<VecR>& pvrefcoords,const VecR& rrefCMass,
                double &rdblMsdMon)
{
    // produces a values for each frame
    int mon;
    int numofmons = pvcoords.size();
    if (numofmons!=pvrefcoords.size()) {
        cout<<"CalcMsdMon::error-> different sizes of coord arrays"<<endl;
        exit(EXIT_FAILURE);
    }
    VecR dr;
    double dummyMSD=0.0;
    for (mon=0;mon<pvcoords.size();mon++){
        dr.x = (pvcoords.at(mon).x - rCMass.x) - (pvrefcoords.at(mon).x - rrefCMass.x);
        dr.y = (pvcoords.at(mon).y - rCMass.y) - (pvrefcoords.at(mon).y - rrefCMass.y);
        dr.z = (pvcoords.at(mon).z - rCMass.z) - (pvrefcoords.at(mon).z - rrefCMass.z);
        dummyMSD += MSQR(dr.x)+MSQR(dr.y)+MSQR(dr.z);
    }
    rdblMsdMon=dummyMSD/(double) numofmons;
}
void CalcEigen(const vector<VecR>& pvcoords,int mysize,VecDoub&d,MatDoub&v){
    // calculates the shape of the molecule fiven the coords pvcoords 
    // and the number of elements that will be consider in the calculation (mysize)
    // Useful in order to calculate the shape of the core and shell monomers
    // Returns two matrices d and v with the EigenValues and eIgenvectors 
    // respectively
    
    int vrows,vcols,dsize;
    const int gyrtensorsize=3;
    
    vrows = v.nrows();
    vcols = v.ncols();
    dsize = d.size();
    
    if ((vrows!=vcols)||(vrows!=dsize)||(vcols!=dsize)||(d.size()!=gyrtensorsize)){
        printf("size of input matrices is wrong!\n");
        exit(2);
    }
    
    MatDoub a(gyrtensorsize,gyrtensorsize);
    
    double  x2,xy,xz,y2,yz,z2;// gyration tensor elements
    int mon;

    x2=0.0; xy=0.0; xz=0.0;
            y2=0.0; yz=0.0;
                    z2=0.0;
    for (mon=0;mon<mysize;mon++){
            x2 += MSQR(pvcoords.at(mon).x - CMass.x);
            xy += (pvcoords.at(mon).x - CMass.x)*(pvcoords.at(mon).y - CMass.y);
            xz += (pvcoords.at(mon).x - CMass.x)*(pvcoords.at(mon).z - CMass.z);

            y2 += MSQR(pvcoords.at(mon).y - CMass.y);
            yz += (pvcoords.at(mon).y - CMass.y)*(pvcoords.at(mon).z - CMass.z);

            z2 += MSQR(pvcoords.at(mon).z - CMass.z);

    }
    x2 /= (double) mysize;
    xy /= (double) mysize;
    xz /= (double) mysize;

    y2 /= (double) mysize;
    yz /= (double) mysize;

    z2 /= (double) mysize;

    a[0][0]=x2; a[0][1]=xy; a[0][2]=xz;
    a[1][0]=xy; a[1][1]=y2; a[1][2]=yz;
    a[2][0]=xz; a[2][1]=yz; a[2][2]=z2;

// Old Jacobi transformation
    //            Jacobi jac(a);
    //            d = jac.d;
    //            v = jac.v;
    //            // order eigenvalues in descenting order
    //            eigsrt(d,&v);
    
    Symmeig sym(a,true);
    d = sym.d;
    v = sym.z;    
    
}
void CalcP2(vector<double>& pvP2,vector< vector<int> >& order,vector<Dend>& myD,int numOfPairs){
    // calculates P2( cos(theta) ) for each per of dendrimers
    static long int P2counter=0;
    int index =0;
	int d1,d2;
	if (P2counter==0) {
		for (int i=0;i<pvP2.size();i++) pvP2.at(i)=0.0;
		P2counter++;
		int i=0;
		for (d1=0;d1<numOfDendrimers-1;d1++){
			for(d2=d1+1;d2<numOfDendrimers;d2++){
				order.at(i).at(0)=d1;
				order.at(i).at(1)=d2;
				i++;
			}
		}
	}
    for (d1=0;d1<numOfDendrimers-1;d1++){
		for(d2=d1+1;d2<numOfDendrimers;d2++){
			pvP2.at(index++)+= MSQR(dot_product(myD.at(d1).laxis,myD.at(d2).laxis));
		}
	}
    
}
int pair(int i,int j,int NN){
	// 2D to 1D index transformation
	return (j-1)+i*(NN-1);	
}
double  MorseForce(MorseParams* pmrs,VecR& dr){
// calculates the force vector between two monomers 
// INPUT:  distance vector dr
// OUTPUT: force vector f 
	double eps,a,b,d;
	double force,rr;
	eps    = pmrs->eps;
	a      = pmrs->a;
	b      = 1.0;
	d      = pmrs->d;
	rr     = sqrt(MSQR(dr.x)+MSQR(dr.y)+MSQR(dr.z));
	if (rr<pmrs->RCut) force  = (2.0*eps*a/b) * exp(-a*(rr-d) )*( exp(-a*(rr-d)) -1 );
	else force=0.0; 

	return force;
}
void FeneForce(VecR& ff,VecR& dr){
// calculates the force vector between two monomers 
// INPUT:  distance vector dr
// OUTPUT: force vector f 
	double K,R,b,L0;
	double force,rr;
	K=0.0;R=0.0;b=0.0;L0=0.0;
	rr =sqrt(MSQR(dr.x)+MSQR(dr.y)+MSQR(dr.z));
	force  = (2.0*K/b)*(rr-L0)/(1-MSQR((rr-L0)/R));
	ff.x   = force*dr.x;
	ff.y   = force*dr.y;
	ff.z   = force*dr.z;
}
double CalcForce(Dend& d1, Dend& d2){
//  Calculates the total force acting from  molecule d2
//  to molecule d1 . The force can be found by integrating the
//  Morse potential energy (MorceForce () )
//  The force is then projected to the line connecting the centre of masses
	int mon1,mon2;
	int typem1,typem2;
	double force;
	double 	result=0.0;
	VecR dr;
	//MorseParams MrsCC,MrsCS,MrsSS;

	MorseParams* pmrs;
//  D7
//  MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.3; MrsCC.RCut = 2.4;
//  MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 2.4;	
//  MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 2.4;

//  D2
//MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.3; MrsCC.RCut = 2.4;
//MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 2.4;	
//MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 2.4;
//	MrsCC.interType = CC ;MrsCC.eps = 0.014;MrsCC.a = 19.2;MrsCC.d = 1.25; MrsCC.RLow = 0.8; MrsCC.RCut = 2.4;	
//	MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 2.4;	
//	MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 2.4;

	for (mon1=0;mon1<d1.moncoords.size();mon1++){
			d1.monforce.at(mon1).x=0.0;
			d1.monforce.at(mon1).y=0.0;
			d1.monforce.at(mon1).z=0.0;
	}
	for (mon2=0;mon2<d2.moncoords.size();mon2++){
			d2.monforce.at(mon2).x=0.0;
			d2.monforce.at(mon2).y=0.0;
			d2.monforce.at(mon2).z=0.0;
	}

	for (mon1=0;mon1<d1.moncoords.size();mon1++){
		typem1 = Montype(numofCoreMons,numofShellMons,d1.id,mon1);
		for (mon2=0;mon2<d2.moncoords.size();mon2++){
			typem2 = Montype(numofCoreMons,numofShellMons,d2.id,mon2);
			
			if (typem1 == typem2){
				if (typem1 == B || typem2 == C) pmrs = &MrsCC;
				else if (typem1 == S)           pmrs = &MrsSS;
			}else if (typem1 != typem2){
				if (typem1 == S || typem2 == S) pmrs = &MrsCS;
				else                            pmrs = &MrsCC;
			}
			if (typem1==-1||typem2==-1) printf("wrong monomer type\n");
			dr.x = d1.moncoords.at(mon1).x - d2.moncoords.at(mon2).x; 
			dr.y = d1.moncoords.at(mon1).y - d2.moncoords.at(mon2).y; 
			dr.z = d1.moncoords.at(mon1).z - d2.moncoords.at(mon2).z; 
			double dumnorma = sqrt(dot_product(dr,dr));
			force=MorseForce(pmrs,dr); // force on mon1 from mon2
			result+=force;
			d1.monforce.at(mon1).x+=force*dr.x/dumnorma;
			d1.monforce.at(mon1).y+=force*dr.y/dumnorma;
			d1.monforce.at(mon1).z+=force*dr.z/dumnorma;

			d2.monforce.at(mon2).x-=force*dr.x/dumnorma;
			d2.monforce.at(mon2).y-=force*dr.y/dumnorma;
			d2.monforce.at(mon2).z-=force*dr.z/dumnorma;
		}
	}
	return result;	
}
double testTotalForce(vector<Dend> & D){
	VecR tf={0.0,0.0,0.0};
	double result=0.0;
	for (int dendid=0;dendid<numOfDendrimers;dendid++){
		for (int mon=0;mon<D.at(dendid).monforce.size();mon++){
			tf.x+=D.at(dendid).monforce.at(mon).x;
			tf.y+=D.at(dendid).monforce.at(mon).y;
			tf.z+=D.at(dendid).monforce.at(mon).z;
		}	
	}
	result+=sqrt(MSQR(tf.x)+MSQR(tf.y)+MSQR(tf.z));
	return result;
}

int Montype(vector<int> & pvcormon,vector<int>& pvshellmon,int dendid,int monid){
	int result=-1;
	if (monid<pvcormon.at(dendid)) result = C;
	else                           result = S;
	return result;
} 

bool OpenCppFileExists(const string& filename)
{
    fstream fin;
    //this will fail if more capabilities to read the 
    //contents of the file is required (e.g. \private\...)
    fin.open(filename.c_str(), ios::in);

    if(fin.is_open())
    {
    fin.close();
    return true;
    }
    fin.close();

    return false;
}
void PotentialSetup(int dt){
	strcpy(fnCC.type,"CC"); strcpy(fnCS.type,"CS"); strcpy(fnBB.type,"BB");
	strcpy(MrsCC.type,"CC");strcpy(MrsCS.type,"CS");strcpy(MrsSS.type,"SS");
    if (dt==1){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.3; MrsCC.RCut = 3.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 3.0;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 3.0;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.7500; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 2.8125; fnBB.R = 0.5625;
    }
    if (dt==2){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.3; MrsCC.RCut = 3.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 3.0;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 3.0;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.7500; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.875; fnBB.R = 0.3750;
    }
    if (dt==3){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.45; MrsCC.RCut = 3.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.50; MrsCS.RLow = 1.25; MrsCS.RCut = 3.0;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 2.00; MrsSS.RLow = 1.70; MrsSS.RCut = 3.0;
        fnCC.interType  = CC ;fnCC.K = 80.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.75; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 80.0; fnBB.L0 = 2.8125; fnBB.R = 0.5625;
    }
    if (dt==4){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.45; MrsCC.RCut = 3.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.50; MrsCS.RLow = 1.25; MrsCS.RCut = 3.0;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 2.00; MrsSS.RLow = 1.70; MrsSS.RCut = 3.0;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.75; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 2.8125; fnBB.R = 0.5625;
    }
    if (dt==5){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 4.8;MrsCC.d = 1.00; MrsCC.RLow = 0.25; MrsCC.RCut = 3.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.50; MrsCS.RLow = 1.2; MrsCS.RCut = 3.0;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 2.00; MrsSS.RLow = 1.7; MrsSS.RCut = 3.0;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.7500; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 2.8125; fnBB.R = 0.5625;
    }
    if (dt==6){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.45; MrsCC.RCut = 3.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.50; MrsCS.RLow = 1.25; MrsCS.RCut = 3.0;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 2.00; MrsSS.RLow = 1.70; MrsSS.RCut = 3.0;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 3.7500; fnCS.R = 0.7500;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.8750; fnBB.R = 0.375;
    }

    if (dt==7){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.0; MrsCC.RCut = 3.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.0; MrsCS.RCut = 3.0;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 0.0; MrsSS.RCut = 3.0;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 2.8125; fnCS.R = 0.5625;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.8750; fnBB.R = 0.3750;
    }
	if (dt==11){
        MrsCC.interType = CC ;MrsCC.eps = 0.014;MrsCC.a = 19.2;MrsCC.d = 1.02; MrsCC.RLow = 0.3; MrsCC.RCut = 2.4;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 2.4;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 2.4;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 2.8125; fnCS.R = 0.5625;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.8750; fnBB.R = 0.3750;
	}
	if (dt==12){
        MrsCC.interType = CC ;MrsCC.eps = 0.714;MrsCC.a = 6.40;MrsCC.d = 1.00; MrsCC.RLow = 0.3; MrsCC.RCut = 1.0;
        MrsCS.interType = CS ;MrsCS.eps = 0.014;MrsCS.a = 19.2;MrsCS.d = 1.25; MrsCS.RLow = 0.8; MrsCS.RCut = 2.4;
        MrsSS.interType = SS ;MrsSS.eps = 0.014;MrsSS.a = 19.2;MrsSS.d = 1.50; MrsSS.RLow = 1.0; MrsSS.RCut = 2.4;
        fnCC.interType  = CC ;fnCC.K = 40.0; fnCC.L0 = 1.8750; fnCC.R = 0.3750;
        fnCS.interType  = CS ;fnCS.K = 20.0; fnCS.L0 = 2.8125; fnCS.R = 0.5625;
        fnBB.interType  = BB ;fnBB.K = 40.0; fnBB.L0 = 1.8750; fnBB.R = 0.3750;
	}
}
