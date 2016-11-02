#include <ctype.h>
#include "system.h"
#include "ran250.h"
#include  <sys/stat.h>
static double RandomLimits(double , double);
void DendSetup() {
    //  numOfMonomersInDend & numOfDendrimers defined in SetParams
    for (int dend = 0; dend < numOfDendrimers; dend++) {
        dendrimer[dend].id = dend;
        dendrimer[dend].numOfMonomers = numOfMonomersInDend;
        dendrimer[dend].gen = G;
        dendrimer[dend].func = F;
        // reset CMass
        dendrimer[dend].CMass.x = 0.0;
        dendrimer[dend].CMass.y = 0.0;
        dendrimer[dend].CMass.z = 0.0;
        dendrimer[dend].Rg.val = 0.0;
        dendrimer[dend].Rg.count = 0;
        dendrimer[dend].Rg.sum = 0.0;
        dendrimer[dend].Rg2.val = 0.0;
        dendrimer[dend].Rg2.count = 0;
        dendrimer[dend].Rg2.sum = 0.0;
        for (int ngen = 0; ngen <= dendrimer[dend].gen; ngen++) {
            dendrimer[dend].numOfMonGen[ngen] =
                    2 * (int) pow((dendrimer[dend].func - 1), ngen);
        }
    }
}
int DendReadCoords(const char inFileName[],int pbc){
	// Read initial configuration file in xyz format
	// UNWRAPPED COORDS
	// reads 1- numOfDendrimers -> monomer.pos	
    // sprintf(inFileName, "%s%d.xyz", fileIniConf, filenum);
    // if pbc  = 0 -> read ini config file with unwrapped positions (read to npos) 
    // if pbc  = 1 -> read ini config file with   wrapped positions (read to  pos)
    
    FILE *pInFile;
    //char *dummy;
    double x, y, z;
    char line[300];
	int dendid,mon;
	int totalmons=0;
	int totaldends=0;
	Dendrimer* pd;

    if ((pInFile = fopen(inFileName, "r")) != NULL) {
    	for (dendid = 0; dendid < numOfDendrimers; dendid++) {
			pd = dendrimer + dendid;
			for (mon = 0; mon < pd->numOfMonomers; mon++) {
				//dummy = fgets(line, 300, pInFile);
				sscanf(line, "%lf%lf%lf", &x, &y, &z);
				if (pbc == 0 ){
					pd->monomer[mon].npos.x = x;
					pd->monomer[mon].npos.y = y;
					pd->monomer[mon].npos.z = z;
				}
				else if (pbc == 1 ){
					pd->monomer[mon].pos.x = x;
					pd->monomer[mon].pos.y = y;
					pd->monomer[mon].pos.z = z;
				}
			
				++totalmons;
			}
			if      (pbc ==0) DendSetpos(pd);
			else if (pbc ==1) DendSetnpos(pd);
			++totaldends;
		}
    }else {
        fprintf(stdout,"DendReadCoords()::inputfile:(%s) not found\n", inFileName);
        exit(3);
    }
	if (totalmons == numOfTotalMonomers){
		fprintf(stdout,"DendReadCoords()::inputfile:(%s) complete read\n", inFileName); 	
		fprintf(stdout,"     total number of dendrimers read: %d\n", totaldends);
		fprintf(stdout,"     total number of monomers   read: %d\n", totalmons);  
	}else{
		fprintf(stdout,"DendReadCoords()::inputfile:(%s) wrong number of monomers read \n", inFileName); 	
		fprintf(stdout,"     total number of dendrimers read: %d\n", totaldends);
		fprintf(stdout,"     total number of monomers   read: %d\n", totalmons);  
		exit(3);
	}
	DendWriteXYZCommon(stdout,0);
	DendWriteXYZCommon(stdout,1);
    fclose(pInFile);

    return 1;
}
int DendSaveConfig(const char outFileName[],int pbc){
	FILE *pOutFile;
	int dendid,mon;
	Dendrimer* pd;
	Monomer *pm;
    if ((pOutFile = fopen(outFileName, "w")) != NULL) {
    	for (dendid = 0; dendid < numOfDendrimers; dendid++) {
        	pd = dendrimer+dendid;
        	for (mon = 0; mon < pd->numOfMonomers; mon++) {
        	    pm = &pd->monomer[mon];
				if (pbc==0)
        	    fprintf(pOutFile, "%16.7E %16.7E %16.7E\n",
        	            (pm->npos).x, (pm->npos).y, (pm->npos).z);
				if (pbc==1)
        	    fprintf(pOutFile, "%16.7E %16.7E %16.7E\n",
        	            (pm->pos).x, (pm->pos).y, (pm->pos).z);
        	}
		}
	}
	fprintf(pOutFile,"energy:%lf\n",RunningEnergy);
	fclose(pOutFile);
	return 1;
}
void DendDefineBondsTypes(Dendrimer *pD) {
    // defines monomer's : id, type, inDend, gen
    static int countmon = 0;
    int first;
    int numOfBonds;
    int Ngf, Ngl, mnext, nbond;

    if (numOfDendrimers > 0) {
        // current dendrimer's (ndend) first monomer
        // bond F-1 is always defines as the bond with the brevious generation
        // exception: central bonded monomers
        first = 0;
        pD->monomer[first].id = first;
        pD->monomer[first].uniqid = countmon++;
        pD->monomer[first].type = B;
        pD->monomer[first].inDend = pD->id;
        pD->monomer[first].gen = 0;
        pD->monomer[first].numOfBonds = pD->func;
        pD->monomer[first].bond[pD->func - 1] = &pD->monomer[first + 1];
        pD->monomer[first].bondInterType[pD->monomer[first].numOfBonds - 1] = BB;

        pD->monomer[first + 1].id = first + 1;
        pD->monomer[first + 1].uniqid = countmon++;
        pD->monomer[first + 1].type = B;
        pD->monomer[first + 1].inDend = pD->id;
        pD->monomer[first + 1].gen = 0;
        pD->monomer[first + 1].numOfBonds = pD->func;
        pD->monomer[first + 1].bond[pD->monomer[first + 1].numOfBonds - 1] = &pD->monomer[first];
        pD->monomer[first + 1].bondInterType[pD->monomer[first + 1].numOfBonds - 1] = BB;

        Ngf = first;
        mnext = first + pD->numOfMonGen[0];
		numOfTotalBonds=0;

        for (int ngen = 0; ngen <= G; ngen++) {
            // Nfg: first monomer's index of generation ngen
            // Ngl-1: last  monomer's index of generation ngen
            // e.g for F = 3
            // ngen=0: Ngf=0;Ngl-1=1;
            // ngen=1: Ngf=2;Ngl-1=5;
            // ngen=2: Ngf=6;Ngl-1=13;
            Ngl = Ngf + pD->numOfMonGen[ngen];
            for (int mon = Ngf; mon < Ngl; mon++) {
                numOfBonds = pD->monomer[mon].numOfBonds;
                for (nbond = 0; nbond < (F - 1); nbond++) {
                    if (ngen < G) {
                        pD->monomer[mnext].id = mnext;
                        pD->monomer[mnext].uniqid = countmon++;
                        pD->monomer[mnext].inDend = pD->id;
                        pD->monomer[mnext].gen = ngen + 1;
                        pD->monomer[mnext].numOfBonds = pD->func;
                        pD->monomer[mnext].bond[numOfBonds - 1] = &pD->monomer[mon];
                        pD->monomer[mon].bond[nbond] = &pD->monomer[mnext];
                        // if ngen == G-1, generation G is build => mnext
                        // belongs to generation G
                        if (ngen == (G - 1)) {
                            pD->monomer[mon].bondInterType[nbond] = CS;
                            pD->monomer[mnext].type = S;
                            pD->monomer[mnext].bondInterType[numOfBonds - 1] = CS;
                        } else {
                            pD->monomer[mnext].type = C;
                            pD->monomer[mon].bondInterType[nbond] = CC;
                            pD->monomer[mnext].bondInterType[numOfBonds - 1] = CC;
                        }
						numOfTotalBonds++;
                    } else if (ngen == G) {
                        pD->monomer[mon].bond[nbond] = NULL;
                        pD->monomer[mon].bondInterType[nbond] = -22;
                    }
                    /*
                                        printf("monomer intertype:%3d\n",pD->monomer[mon].bondInterType[nbond]);
                     */
                    ++mnext;
                }
            }
            Ngf = Ngl;
        }
    }
}

void DendCalcCMass (Dendrimer *pd) {
	// USES NPOS
    // unfold function has to be called before calculating CMass;
    int mon;
    // Centre of Mass Calculation
    pd->CMass.x = 0.0;
    pd->CMass.y = 0.0;
    pd->CMass.z = 0.0;

    for (mon = 0; mon < pd->numOfMonomers; mon++) {
        pd->CMass.x += pd->monomer[mon].npos.x;
        pd->CMass.y += pd->monomer[mon].npos.y;
        pd->CMass.z += pd->monomer[mon].npos.z;
    }
    pd->CMass.x /= (double) pd->numOfMonomers;
    pd->CMass.y /= (double) pd->numOfMonomers;
    pd->CMass.z /= (double) pd->numOfMonomers;
}

void DendCalcRg(Dendrimer * pD) {
    // unfold function has to be called before calculating CMass;
    static long int rgcount = 0;
    int mon;
    double rg, rg2;

    rg2 = 0.0;

    for (mon = 0; mon < pD->numOfMonomers; mon++) {
        rg2 += SQR(pD->monomer[mon].npos.x - pD->CMass.x) +
               SQR(pD->monomer[mon].npos.y - pD->CMass.y) +
               SQR(pD->monomer[mon].npos.z - pD->CMass.z);
    }

    rg2 /= (double) pD-> numOfMonomers;
    rg = sqrt(rg2);

    pD->Rg2.val = rg2;
    pD->Rg.val = rg;

    if (rgcount == 0) {
        pD->Rg.sum = 0.0;
        pD->Rg.count = 0;
        pD->Rg2.sum = 0.0;
        pD->Rg2.count = 0;
    }

    pD->Rg2.sum += rg2;
    pD->Rg2.count++;
    pD->Rg.sum += rg;
    pD->Rg.count++;
    rgcount++;

}

double *DendGyrTensor(Dendrimer *pD){

	double *gyr;
	const int gyrsize=3;	

	AllocMat(gyr,gyrsize,gyrsize,double);
	int mon;

	double 	x2=0.0,xy=0.0,xz=0.0,y2=0.0,yz=0.0,z2=0.0;
	int numOfMonomersInDend = pD->numOfMonomers;

	for (mon=0;mon<pD->numOfMonomers;mon++){

		x2 += SQR((pD->monomer)[mon].npos.x - pD->CMass.x);
		xy +=    ((pD->monomer)[mon].npos.x - pD->CMass.x)*((pD->monomer)[mon].npos.y - pD->CMass.y);
		xz +=    ((pD->monomer)[mon].npos.x - pD->CMass.x)*((pD->monomer)[mon].npos.z - pD->CMass.z);
		
		y2 += SQR((pD->monomer)[mon].npos.y - pD->CMass.y);
		yz +=    ((pD->monomer)[mon].npos.y - pD->CMass.y)*((pD->monomer)[mon].npos.z - pD->CMass.z);
	
		z2 += SQR((pD->monomer)[mon].npos.z - pD->CMass.z);
	
	}
	
	x2 /= (double) numOfMonomersInDend;
	xy /= (double) numOfMonomersInDend;
	xz /= (double) numOfMonomersInDend;
	
	y2 /= (double) numOfMonomersInDend;
	yz /= (double) numOfMonomersInDend;
	
	z2 /= (double) numOfMonomersInDend;
	// gyration tensor matrix
	// A(i,j)= gyr[i*cols+j)
	gyr[0]=x2; gyr[1]=xy; gyr[2]=xz;
	gyr[3]=xy; gyr[4]=y2; gyr[5]=yz;
	gyr[6]=xz; gyr[7]=yz; gyr[8]=z2;
#ifdef CALC_EIGEN
	
#endif
	return gyr;
}

double DendMaxRad(const Dendrimer *pd) {
    return fnBB.L0 + 2 * (pd->gen - 1) * fnCC.L0 + 2 * fnCS.L0;
}
void DendWriteProps(const Dendrimer *pD, FILE * pOutFile) {
    int monBondInterType;
    fprintf(pOutFile, "Properties of Dendrimer(%-d)\n", pD->id);
    fprintf(pOutFile, "*****************************\n");
    fprintf(pOutFile, "Total Number of Monomers:  %3d\n",
            pD->numOfMonomers);
    fprintf(pOutFile, "id:%-4d gen:%-3d func:%-3d\n", pD->id,
            pD->gen, pD->func);
    fprintf(pOutFile, "num Of mon in each generation:\n");
    for (int gen = 0; gen <= pD->gen; gen++) {
        fprintf(pOutFile, "G%3d :%3d\t", gen, pD->numOfMonGen[gen]);
    }
    fprintf(pOutFile, "\n");
    fprintf(pOutFile, "     Monomer Properties\n");
    fprintf(pOutFile, " mon    bonded with\t\tposition(x,  y,  z)\n");
    for (int mon = 0; mon < pD->numOfMonomers; mon++) {
        fprintf(pOutFile, "  %-4d ", pD->monomer[mon].id);
        for (int nbond = 0; nbond < pD->monomer[mon].numOfBonds; nbond++) {
            monBondInterType = pD->monomer[mon].bondInterType[nbond];
            if (pD->monomer[mon].bond[nbond] != NULL) {
                fprintf(pOutFile, "%3d", pD->monomer[mon].bond[nbond]->id);
            } else {
                fprintf(pOutFile, "%3d", -2);
            }
            fprintf(pOutFile, "%2s\t",
                    ((monBondInterType == 0) ? "(CC)" :
                    ((monBondInterType == 1) ? "(CS)" :
                    ((monBondInterType == 2) ? "(BB)" : "X"))));
        }
        fprintf(pOutFile, "%13.6e\t%13.6e\t%13.6e",
//                pD->monomer[mon].npos.x, pD->monomer[mon].npos.y, pD->monomer[mon].npos.z);
                pD->monomer[mon].pos.x, pD->monomer[mon].pos.y, pD->monomer[mon].pos.z);
        fprintf(pOutFile, "\n");
    }
    fprintf(pOutFile, "-------------------------------\n");
}
void DendWriteXYZCommon(FILE *pOutFile,int pbc) {
	static long int framenum=1;
    int dendid, mon;
    int monmark = 0;
    int monspc;
    Dendrimer *pd;
    Monomer *pm;
	//if (framenum>1) fprintf(pOutFile,"%-ld\n",framenum);
    for (dendid = 0; dendid < numOfDendrimers; dendid++) {
        pd = dendrimer+dendid;
        for (mon = 0; mon < pd->numOfMonomers; mon++) {
            pm = &pd->monomer[mon];
            monspc = pm->type + 1;
			if (pbc==0)
            fprintf(pOutFile, "%16.7E %16.7E %16.7E %4.1f %1d %1d\n",
                    (pm->npos).x, (pm->npos).y, (pm->npos).z, 0.5, monspc,monmark);
			if (pbc==1)
            fprintf(pOutFile, "%16.7E %16.7E %16.7E %4.1f %1d %1d\n",
                    (pm->pos).x, (pm->pos).y, (pm->pos).z, 0.5, monspc,monmark);
        }
    }
    fprintf(pOutFile, "\n\n");
	framenum++;
}
void DendWriteCMassXYZ(FILE * pOutFile,int pbc){
	static long int framenum=1;
    int dendid;
    Dendrimer *pd;
	fprintf(pOutFile,"%-5d\n%ld\n",numOfDendrimers,framenum);
    for (dendid = 0; dendid < numOfDendrimers; dendid++) {
        pd = dendrimer+dendid;
		fprintf(pOutFile, "%5d\t%16.7E\t%16.7E\t%16.7E\n",
                    dendid,(pd->CMass).x,(pd->CMass).y,(pd->CMass).z);
    }
	framenum++;
	
}
void DendWriteVTKframe(){
	static long int framenum=1;
	FILE *pFile;
	char filename[120];
	Dendrimer *pd;
	VecR cm;
	sprintf(filename,"gyrtensor%06ld.vtk",framenum);
	pFile=fopen(filename,"w");
	
	fprintf(pFile,"# vtk DataFile Version 3.0\n");
	fprintf(pFile,"Random data to test tensors\n");
	fprintf(pFile,"ASCII\n");
	fprintf(pFile,"DATASET UNSTRUCTURED_GRID\n");
	fprintf(pFile,"POINTS %d float\n",numOfDendrimers);
	for (int dendid=0;dendid<numOfDendrimers;dendid++){
		pd = dendrimer + dendid;
		VCopy(cm,(pd->CMass));
		PBCAll(cm);
		
	}
}
void WriteBondVmd (const char outfilename[]){
	FILE *pFile;
    if ((pFile = fopen(outfilename, "w")) != NULL) {
	    fprintf(pFile,"$sel setbonds {\n");
		for (int dendid=0;dendid<numOfDendrimers;dendid++){
			Dendrimer *pd=dendrimer+dendid;
			for (int mon=0;mon<pd->numOfMonomers;mon++){
				Monomer *pm=(pd->monomer)+mon;
	    	    fprintf(pFile,"{ %5d ",pm->uniqid);
	    	    for (int nbond=0;nbond<pm->numOfBonds;nbond++){
	    	        if (pm->bond[nbond]!=NULL){
	    	            fprintf(pFile," %5d ",(pm->bond[nbond])->uniqid);
	    	        }
	    	    }
	    	    fprintf(pFile,"}\n");
	    	}
		}
	    fprintf(pFile,"}\n");
	}
	else{
		return;
	}
	fclose(pFile);

}
void DendWriteinlmp(const char outfilename[],const char outfilename2[],const int pbc){
	FILE *pOutFile;
	Dendrimer *pd;
	Monomer   *pm,*pb;
	int countmol,countatoms,countbonds,atomtype,bondtype;
	double mass=1.0;
	countatoms = 0 ;
	countbonds = 0 ;
	countmol  = 0;
    if ((pOutFile = fopen(outfilename, "w")) != NULL) {
		fprintf(pOutFile,"#Topology file--Num dendrimers:%d, pottype:%d,Total energy:%lf\n\n",numOfDendrimers,PotType,EnergySystem());
		fprintf(pOutFile,"\t%d\tatoms\n\t%d\tbonds\n\t%d\tdihedrals\n\t%d\timpropers\n\n",numOfTotalMonomers,numOfDendrimers*(numOfTotalBonds+1),0,0);
		fprintf(pOutFile,"\t%d\t atom types\n\t%d\t bond types\n\n",NTYPE,N_TYPE_FENE);
		fprintf(pOutFile,"\t%8.3f\t%8.3f\txlo xhi\n\t%8.3f\t%8.3f\tylo yhi\n\t%8.3f\t%8.3f\tzlo zhi\n\n\n",
							-Box.x/2.0,Box.x/2.0,-Box.y/2.0,Box.y/2.0,-Box.z/2.0,Box.z/2.0);
		fprintf(pOutFile,"Masses\n\n");
		for (int i=0;i<NTYPE;i++){
			fprintf(pOutFile,"\t%d\t%4.1f\n",i+1,mass);
		}
		fprintf(pOutFile,"\nAtoms\n\n");

//	    ATOM COORDINATES
		for (int dendid = 0; dendid < numOfDendrimers; dendid++) {
			countmol++;
    	    pd = dendrimer+dendid;
    	    for (int mon = 0; mon < pd->numOfMonomers; mon++) {
				countatoms++;
    	        pm = &pd->monomer[mon];
				if ((pm->type == C)||(pm->type == B)) atomtype = C ;
				else if (pm->type == S)               atomtype = S ;
				fprintf(pOutFile,"\t\t%d\t%d\t%d\t\t",countatoms,countmol,atomtype);
				if (pbc==0)
    	        fprintf(pOutFile, "%17.8E\t%17.8E\t%17.8E\t\t%1d %1d %1d\n",
    	                (pm->npos).x, (pm->npos).y, (pm->npos).z, 0, 0, 0);
				else if (pbc==1)
    	        fprintf(pOutFile, "%17.8E\t%17.8E\t%17.8E\t\t%1d %1d %1d\n",
    	                (pm->pos).x, (pm->pos).y, (pm->pos).z, 0, 0, 0);
    	    }
    	}
//		BONDS		
		fprintf(pOutFile,"\n\nBonds\n\n");
		for (int dendid = 0; dendid < numOfDendrimers; dendid++) {
    	    pd = dendrimer+dendid;
    	    for (int mon = 0; mon < pd->numOfMonomers; mon++) {
    	        pm = &pd->monomer[mon];
    	        for (int nbond = 0; nbond < pm->numOfBonds; nbond++) {
    	            pb       = pm->bond[nbond];
					bondtype = (pm->bondInterType)[nbond]+1;
    	            if (pb != NULL) {
    	                if (pm->id < pb->id) {
							countbonds++;
    	                    fprintf(pOutFile, "\t%d\t%d\t%d\t%d\t\n",
    	                            countbonds,bondtype,dendid*(pd->numOfMonomers)+pm->id+1 ,dendid*(pd->numOfMonomers)+pb->id +1);
    	                }
    	            }
    	        }
    	    }
    	}
		fprintf(pOutFile,"\n");
		fclose(pOutFile);
    	if ((pOutFile = fopen(outfilename2, "w")) != NULL) {
			fprintf(pOutFile,"pair_style morse %lf\n",3.0);
			fprintf(pOutFile,"pair_coeff %d %d %9.4f %9.4f %9.4f %9.4f\n",C,C,MrsCC.eps,MrsCC.a,MrsCC.d,MrsCC.RCut);
			fprintf(pOutFile,"pair_coeff %d %d %9.4f %9.4f %9.4f %9.4f\n",C,S,MrsCS.eps,MrsCS.a,MrsCS.d,MrsCS.RCut);
			fprintf(pOutFile,"pair_coeff %d %d %9.4f %9.4f %9.4f %9.4f\n\n",S,S,MrsSS.eps,MrsSS.a,MrsSS.d,MrsSS.RCut);
			fprintf(pOutFile,"pair_modify shift yes\n");
			fprintf(pOutFile,"bond_style fene/expand\n");
			fprintf(pOutFile,"special_bonds fene\n\n");
			fprintf(pOutFile,"bond_coeff %d %9.4f %9.4f %9.4f %9.4f %9.4f\n",CC+1,2.0*fnCC.K,fnCC.R,0.0,0.0,fnCC.L0);
			fprintf(pOutFile,"bond_coeff %d %9.4f %9.4f %9.4f %9.4f %9.4f\n",CS+1,2.0*fnCS.K,fnCS.R,0.0,0.0,fnCS.L0);
			fprintf(pOutFile,"bond_coeff %d %9.4f %9.4f %9.4f %9.4f %9.4f\n\n",BB+1,2.0*fnBB.K,fnBB.R,0.0,0.0,fnBB.L0);
		}
	}
}
void DendWriteBonds(FILE* pFile) {
    int countmon;
    int dendid;
    Dendrimer *pd;
    int mon, nbond;
    double bondrad = 0.2;
    Monomer *pm, *pb;

    int bonddraw = 1;

    countmon = 0;
    for (dendid = 0; dendid < numOfDendrimers; dendid++) {
        pd = &dendrimer[dendid];
        countmon = (pd->id) * (pd->numOfMonomers);
        for (mon = 0; mon < pd->numOfMonomers; mon++) {
            pm = &pd->monomer[mon];
            for (nbond = 0; nbond < pm->numOfBonds; nbond++) {
                pb = pm->bond[nbond];
                if (pb != NULL) {
                    if (pm->id < pb->id) {
                        fprintf(pFile, "%3d %3d  %4.2f  %1d\n",
                                pm->id + countmon, pb->id + countmon,
                                bondrad, bonddraw);
                    }
                }
            }
        }
    }

}
void DendInitConfig(int dendid, const VecR *pvecr) {
    // Creates RANDOM configuration of dendrimer
    // The first(central) monomer is positioned at
    // pvecr
    //printf("Building initial configurations for Dendrimer:%-d\n", dendid);
    int first, numOfBonds;
    int mnext, nbond, ngen, Ngf, Ngl;
    double radius;
    int counter;
    VecR vecr;

    nbond = 0;
    ngen = 0;
    Dendrimer *pD;

    pD = dendrimer+dendid;
    //
    //place gen=0 monomers, central monomers
    radius= RandomLimits(fnBB.L0 - fnBB.R, fnBB.L0 + fnBB.R);//radius = fnBB.L0;
    // set the position of the first monomer to
    // the centre of the Box
    first = 0;
    pD->monomer[first].pos.x = pvecr->x; // (double) (Box.x/2);
    pD->monomer[first].pos.y = pvecr->y; // (double) (Box.y/2);
    pD->monomer[first].pos.z = pvecr->z; //(double) (Box.z/2);
    RandomVector(&vecr, radius);
    pD->monomer[first + 1].pos.x = pD->monomer[first].pos.x + vecr.x;
    pD->monomer[first + 1].pos.y = pD->monomer[first].pos.y + vecr.y;
    pD->monomer[first + 1].pos.z = pD->monomer[first].pos.z + vecr.z;

    Ngf = first;
    mnext = first + dendrimer[dendid].numOfMonGen[0];

    for (ngen = 0; ngen < G; ngen++) {
        Ngl = Ngf + dendrimer[dendid].numOfMonGen[ngen];
        for (int mon = Ngf; mon < Ngl; mon++) {
            if ((pD->monomer[mon].type == B || pD->monomer[mon].type == C)&& ngen == (G - 1)){
                radius= RandomLimits(fnCS.L0 - fnCS.R, fnCS.L0 + fnCS.R);//radius = fnCS.L0;   
            }else{
                radius= RandomLimits(fnCC.L0 - fnCC.R, fnCC.L0 + fnCC.R);//radius = fnCC.L0;
            }
            numOfBonds = pD->monomer[mon].numOfBonds;
            for (nbond = 0; nbond < (numOfBonds - 1); nbond++) {
                counter++;
                RandomVector(&vecr, radius);
                pD->monomer[mnext].pos.x = pD->monomer[mon].pos.x + vecr.x;
                pD->monomer[mnext].pos.y = pD->monomer[mon].pos.y + vecr.y;
                pD->monomer[mnext].pos.z = pD->monomer[mon].pos.z + vecr.z;
                ++mnext;
            }
        }
        Ngf = Ngl;
    }
	DendSetnpos(pD);
}

static double RandomLimits(double dblmin,double dblmax){
    // returns a random double number from (dblmin,dblmax)
    return ((dblmax-dblmin)*dr250() + dblmin );
}

int DendCheckBox() {
    int mon;
    int result = FALSE;
    Dendrimer *pd;
    Monomer *pm;
    

    for (int dendid=0; dendid<numOfDendrimers;dendid++){
        pd = &dendrimer[dendid];
        for (mon = 0; mon < pd->numOfMonomers; mon++) {
            pm = &pd->monomer[mon];
            if (pm->pos.x > Box.x / 2.0 || pm->pos.x < -Box.x / 2.0) {
                result = FALSE;
                WriteVecR(pm->pos,stdout);
                break;
            } else if (pm->pos.y > Box.y / 2.0 || pm->pos.y < -Box.y / 2.0) {
                result = FALSE;
                WriteVecR(pm->pos,stdout);
                break;
            } else if (pm->pos.z > Box.z / 2.0 || pm->pos.z < -Box.z / 2.0) {
                result = FALSE;
                WriteVecR(pm->pos,stdout);
                break;
            } else {
                result = TRUE;
            }

        }
    }
    return result;
}


void DendSetnpos(Dendrimer *pd){
    for(int mon=0;mon<pd->numOfMonomers;mon++){
        // allocate pointer to position in case is NULLified
        pd->monomer[mon].npos.x = pd->monomer[mon].pos.x;
        pd->monomer[mon].npos.y = pd->monomer[mon].pos.y;
        pd->monomer[mon].npos.z = pd->monomer[mon].pos.z;
    }
}
void DendSetpos(Dendrimer *pd){
    for(int mon=0;mon<pd->numOfMonomers;mon++){
        pd->monomer[mon].pos.x = pd->monomer[mon].npos.x;
        pd->monomer[mon].pos.y = pd->monomer[mon].npos.y;
        pd->monomer[mon].pos.z = pd->monomer[mon].npos.z;
    }
}
/*
void DendUnfoldXYZ(){
    // unfolds the pbc coordinates after initializing them to
    // the pbc position, reference monomer has to be inside the box
    int ndend;
    int nmon;
    Dendrimer *pd;
    VecR vecr;

    for (ndend=0;ndend<numOfDendrimers;ndend++){
        pd=dendrimer+ndend;

        DendSetnpos(pd);
        vecr.x = pd->monomer[0].pos.x;
        vecr.y = pd->monomer[0].pos.y;
        vecr.z = pd->monomer[0].pos.z;
        
        for (nmon=0;nmon<pd->numOfMonomers;nmon++){
              DendPBCunwrap(&pd->monomer[nmon].npos,&vecr);
        }
    }

}
void DendPBCunwrap(VecR* pos, const VecR *refpos){
    pos->x -= lround( ( pos->x - refpos->x ) / Box.x )*(Box.x);
    pos->y -= lround( ( pos->y - refpos->y ) / Box.y )*(Box.y);
    pos->z -= lround( ( pos->z - refpos->z ) / Box.z )*(Box.z);
}
*/

/*
void DendSetnposAll(){
    // set *npos of all monomers to point to pos
    // no pbc no need to change functions
    for (int dendid=0;dendid<numOfDendrimers;dendid++){
        for (int mon=0;mon<dendrimer[dendid].numOfMonomers;mon++){
            dendrimer[dendid].monomer[mon].npos = &dendrimer[dendid].monomer[mon].pos;
        }
    }
}

void DendUnSetnposAll(){
    // set *npos of all monomers to point to pos
    // no pbc no need to change functions
    for (int dendid=0;dendid<numOfDendrimers;dendid++){
        for (int mon=0;mon<dendrimer[dendid].numOfMonomers;mon++){
            dendrimer[dendid].monomer[mon].npos = NULL;
        }
    }
}
*/
/*
for (;;) {
                    counter++;
                    RandomVector(&vecr, radius);
                    pD->monomer[mnext].pos.x = pD->monomer[mon].pos.x + vecr.x;
                    pD->monomer[mnext].pos.y = pD->monomer[mon].pos.y + vecr.y;
                    pD->monomer[mnext].pos.z = pD->monomer[mon].pos.z + vecr.z;
                    ok = 1;
                    for (i = 0; i < mnext; i++) {
                        //                        printf("i=%-3d\n",i);
                        dr.x = pD->monomer[mnext].pos.x - pD->monomer[i].pos.x;
                        dr.y = pD->monomer[mnext].pos.y - pD->monomer[i].pos.y;
                        dr.z = pD->monomer[mnext].pos.z - pD->monomer[i].pos.z;
                        r2 = SQR(dr.x) + SQR(dr.y) + SQR(dr.z);
                        if (radius > 1.0) tol = 1.0;
                        else tol = SQR(radius);
                        if (r2 < tol) {
                            ok = 1;
                            break;
                        } else{
                            ok = 1;
                        }
                    }
                    if (ok == 1) break;
                }
*/
/*void DendWriteXYZ() {
    static int count = 0;
    char outFileName[130];
    int dendid;
    Dendrimer *pd;
    if (count == 0) {
        for (dendid = 0; dendid < numOfDendrimers; dendid++) {
            pd = &dendrimer[dendid];
            sprintf(outFileName, "MC_d%1dG%1dF%1d.xyz",
                    pd->id, pd->gen, pd->func);
            pOutputFile[dendid] = fopen(outFileName, "w");
            fprintf(pOutputFile[dendid], "%-5d \n \n", pd->numOfMonomers);
        }
        count++;
    }
    for (dendid = 0; dendid < numOfDendrimers; dendid++) {
        pd = &dendrimer[dendid];
        for (int mon = 0; mon < pd->numOfMonomers; mon++) {
            fprintf(pOutputFile[dendid], "G%-1d\t%8.3f\t%8.3f\t%8.3f\n",
                    pd->monomer[mon].gen, pd->monomer[mon].pos.x,
                    pd->monomer[mon].pos.y, pd->monomer[mon].pos.z);
        }
        fprintf(pOutputFile[dendid], "\n\n");
    }

}*/
