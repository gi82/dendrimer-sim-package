#include "system.h"
#include "preprocessor.h"

double InterDend(Dendrimer *pd1, Dendrimer *pd2){
    // CALCULATES THE ENERGY OF ALL MONOMERS OF DENDRIMER 1 (pd1)
    // WITH MONOMERS OF DENDRIMER 2 (pd1)
    // *** ONLY MORSE POTENTIAL ***
    int     mon1, mon2;
    Monomer *pm1, *pm2;
    int     typem1,typem2;
    VecR    dr;
    double  rsq, result;
    result =0.0;

#ifdef LUT_MORSE
	int interType=-1;
#else
    MorseParams* pmrs =NULL;
#endif


    // LOOP OVER ALL MONOMERS OF DENDRIMERS 1 & 2
    for (mon1 =0;mon1<pd1->numOfMonomers;mon1++){
        pm1 = &pd1->monomer[mon1];
        typem1 = pm1->type;
        for (mon2 = 0 ;mon2 < pd2->numOfMonomers;mon2++){
                pm2 = &pd2->monomer[mon2];
                dr.x = pm1->pos.x - pm2->pos.x;
                dr.y = pm1->pos.y - pm2->pos.y;
                dr.z = pm1->pos.z - pm2->pos.z;

                typem2 = pm2->type;

#if USE_PBC==1
                PBCAll(dr);
#endif
                rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
#ifdef LUT_MORSE
				if (typem1 == typem2){
					if (typem1 == B || typem2 == C) {interType = CC;}
					else if (typem1 == S)           {interType = SS;}
				}else if (typem1 != typem2){         
					if (typem1 == S || typem2 == S) {interType = CS;}
					else                            {interType = CC;}
				}
				result += LUTMorseValue (interType,rsq);
#else
				if (typem1 == typem2){
					if (typem1 == B || typem2 == C) {pmrs = &MrsCC;}
					else if (typem1 == S)           {pmrs = &MrsSS;}
				}else if (typem1 != typem2){                       
					if (typem1 == S || typem2 == S) {pmrs = &MrsCS;}
					else                            {pmrs = &MrsCC;}
				}
                if (rsq<SQR(pmrs->RCut)){
                    result += Morse(pmrs,sqrt(rsq))-Morse(pmrs,pmrs->RCut);
                }
#endif
        }
    }
    return result;
}

double IntraDend(Dendrimer *pd){
    // calculate internal energy of dendrimer pd
    int     mon1, mon2, nbond;
    Monomer *pm1, *pm2, *pb;
    double  result, rsq;
    VecR    dr;
    int     typem1, typem2;;
    result =0;

    // BONDS Fene + Morse
#if defined(LUT_MORSE) || defined(LUT_FENE)
	int interType = -1;
#endif
#ifndef LUT_MORSE
    MorseParams* pmrs =NULL;
#endif
#ifndef LUT_FENE
    FeneParams *pfene=NULL;
#endif

    for (mon1 =0; mon1 <pd->numOfMonomers ; mon1++){
        pm1 = &pd->monomer[mon1];
        typem1 = pm1->type;
        for (nbond =0; nbond <pm1->numOfBonds; nbond++){
            pb = pm1->bond[nbond];
            if ((pb != NULL)&&(pb->id>pm1->id)){
                dr.x = pm1->pos.x-pb->pos.x;
                dr.y = pm1->pos.y-pb->pos.y;
                dr.z = pm1->pos.z-pb->pos.z;
#if USE_PBC==1
                PBCAll(dr);
#endif
                rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
#ifdef LUT_FENE
                if      (pm1->bondInterType[nbond]==CC) {interType = CC;}
                else if (pm1->bondInterType[nbond]==CS) {interType = CS;}
                else if (pm1->bondInterType[nbond]==BB) {interType = BB;}
				result += LUTFeneValue (interType, rsq);
#else
                if      (pm1->bondInterType[nbond]==CC) { pfene = &fnCC;}
                else if (pm1->bondInterType[nbond]==CS) { pfene = &fnCS;}
                else if (pm1->bondInterType[nbond]==BB) { pfene = &fnBB;}
                else pfene = NULL;
                result += Fene(pfene,sqrt(rsq));
#endif

            }
        }
        // Morse Interaction
        // for (mon2 =mon1+1 ;mon2 <pd->numOfMonomers;mon2++){}
        // gives the same results in almost the same time
    }
    // Morse
#if ONLY_MORSE==1
    // calculate only morse interaction for testing purposes
    result = 0.0;
#endif
    
    for (mon1 =0; mon1 <pd->numOfMonomers-1; mon1++){
        pm1 = &pd->monomer[mon1];
        typem1 = pm1->type;
        for (mon2 =mon1+1 ;mon2 <pd->numOfMonomers;mon2++){
                pm2 = &pd->monomer[mon2];
                typem2 = pm2->type;
                dr.x = pm1->pos.x - pm2->pos.x;
                dr.y = pm1->pos.y - pm2->pos.y;
                dr.z = pm1->pos.z - pm2->pos.z;
#if USE_PBC==1
                PBCAll(dr);
#endif
                rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
#ifdef LUT_MORSE
				if (typem1 == typem2){
					if (typem1 == B || typem2 == C) { interType = CC;}
					else if (typem1 == S)           { interType = SS;}
				}else if (typem1 != typem2){          
					if (typem1 == S || typem2 == S) { interType = CS;}
					else                            { interType = CC;}
				}
				result += LUTMorseValue (interType,rsq);
#else
				if (typem1 == typem2){
					if (typem1 == B || typem2 == C) { pmrs = &MrsCC; }
					else if (typem1 == S)           { pmrs = &MrsSS; }
				}else if (typem1 != typem2){                         
					if (typem1 == S || typem2 == S) { pmrs = &MrsCS; }
					else                            { pmrs = &MrsCC; }
				}
                if (rsq<SQR(pmrs->RCut)){
                    result += Morse(pmrs,sqrt(rsq))-Morse(pmrs,pmrs->RCut);
                }
#endif
        }
    }
    return result;
}
double PotDendDend(Dendrimer *pd1, Dendrimer *pd2){
    double result =0.0;
    int mon;
    for (mon =0;mon< pd1->numOfMonomers;mon++){
        result += PotMonNextDend(&pd1->monomer[mon],&pd1->monomer[mon].pos,pd2);
    }
    return result;
}
double PotMonNextDend(Monomer *pm, VecR *monpos, Dendrimer *pd){
    // Energy contribution of monomer *pm with dendrimer *pd
    // Only Morse Contributions
    int    typem,  typemnext;
    double rsq, result =0.0;
    VecR   dr;

    Monomer *pmnext;

#ifdef LUT_MORSE
	int interType = -1;
#else
    MorseParams* pmrs = NULL;
#endif


    if (pm!=NULL){
        if (pm->inDend != pd->id){
            typem = pm->type;
            for (int monnext =0; monnext <pd->numOfMonomers; monnext++){
                pmnext = &pd->monomer[monnext];
                dr.x=monpos->x-pmnext->pos.x;
                dr.y=monpos->y-pmnext->pos.y;
                dr.z=monpos->z-pmnext->pos.z;
#if USE_PBC==1
                PBCAll(dr);
#endif
                rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                typemnext = pmnext->type;
#ifdef LUT_MORSE
                if (typem == typemnext){
                    if (typem == B || typemnext == C) { interType = CC;}
                    else if (typem == S)              { interType = SS;}
                }else if (typem != typemnext){
                    if (typem == S || typemnext == S) { interType = CS;}
                    else                              { interType = CC;}
                }
				result += LUTMorseValue (interType,rsq);
#else
                if (typem == typemnext){
                    if (typem == B || typemnext == C) { pmrs = &MrsCC; }
                    else if (typem == S)              { pmrs = &MrsSS; }
                }else if (typem != typemnext){
                    if (typem == S || typemnext == S) { pmrs = &MrsCS; }
                    else                              { pmrs = &MrsCC; }
                }
                if (rsq<SQR(pmrs->RCut)){
                    result += Morse(pmrs,sqrt(rsq))-Morse(pmrs,pmrs->RCut);
                }
#endif
            }
        }else{
            result = 0.0;
        }
    }
    return result;
}

double CelldESingleMonMove(Monomer* pmon1, VecR *monpos){
// Calculates the total Morse Energy difference for the trial monomer move
// of monomer pm from pm->pos to monpos
    
    double   result=0.0;
    int      cindex, nei;
    Cell     *pcell1,*pcell2;
    Monomer  *pmon2;
    int      typem1, typem2;
    VecR     dr;
    double   rsq;
    int      mon2;

#ifdef LUT_MORSE
	int      interType = -1;
#else
    MorseParams *pmrs=NULL;
#endif


    cindex = pmon1->cellIndex;

    pcell1 = cells + cindex;
    typem1 = pmon1->type;
 
    for (nei=0;nei<pcell1->numOfNeigh;nei++){
        pcell2 = cells + (pcell1->neighbours)[nei];
        pmon2  = pcell2->first;
        for (mon2=0;mon2<pcell2->numOfElements;mon2++){
            if ((pmon1!=NULL)&&(pmon2!=NULL)){
                if (pmon1!= pmon2){
                    dr.x = monpos->x - (pmon2->pos).x;
                    dr.y = monpos->y - (pmon2->pos).y;
                    dr.z = monpos->z - (pmon2->pos).z;
#if USE_PBC==1
                    PBCAll(dr);
#endif
                    rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                    typem2=pmon2->type;
#ifdef LUT_MORSE
                    if (typem1 == typem2){
                        if (typem1 == B || typem2 == C) {interType = CC;}
                        else if (typem1 == S)           {interType = SS;}
                    }else if (typem1 != typem2){
                        if (typem1 == S || typem2 == S) {interType = CS;}
                        else                            {interType = CC;}
                    }
					result += LUTMorseValue (interType, rsq);
#else
                    if (typem1 == typem2){
                        if (typem1 == B || typem2 == C) { pmrs = &MrsCC;}
                        else if (typem1 == S)           { pmrs = &MrsSS;}
                    }else if (typem1 != typem2){
                        if (typem1 == S || typem2 == S) { pmrs = &MrsCS;}
                        else                            { pmrs = &MrsCC;}
                    }
                	if (rsq<SQR(pmrs->RCut)){
                	    result += Morse(pmrs,sqrt(rsq))-Morse(pmrs,pmrs->RCut);
                	}
#endif
                }                
            }
            
            pmon2 = pmon2->next;
        }
    }
    return result;
}

double  PotMonIntraDend(Monomer* pm, VecR *monpos,int mb){
    // calculates total interaction between monomer [mon] and the rest of the
    // monomer of the SAME dendrimer [dend]
    // VecR trialpos is used
    int    dend, nbond;
    int    typem,monnext,typemnext;
    double rsq;    // squared distance
    Dendrimer *pd;
    Monomer   *pb, *pmnext; // current and bonded monomer
    VecR       dr;

#if defined(LUT_MORSE) || defined(LUT_FENE)
	int interType      = -1;
#endif
#ifndef LUT_MORSE
    MorseParams *pmrs  = NULL;
#endif
#ifndef LUT_FENE
    FeneParams  *pfene = NULL;
#endif
    
    dend = pm->inDend;
    pd   = dendrimer+dend;

    // BONDED MONOMERS - BONDED MONOMERS
    double fene = 0.0;
    for (nbond = 0; nbond < pm->numOfBonds ; nbond++){
        pb = pm->bond[nbond];
        if (pb != NULL){//
            dr.x = monpos->x-pb->pos.x;
            dr.y = monpos->y-pb->pos.y;
            dr.z = monpos->z-pb->pos.z;
#if USE_PBC==1
            PBCAll(dr);
#endif
            rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
#ifdef LUT_FENE
			if      (pm->bondInterType[nbond]==CC) { interType = CC;}
			else if (pm->bondInterType[nbond]==CS) { interType = CS;}
			else if (pm->bondInterType[nbond]==BB) { interType = BB;}
			else {
				printf("Wrong Fene\n");
			}
			fene += LUTFeneValue (interType, rsq);
#else
			if      (pm->bondInterType[nbond]==CC) { pfene = &fnCC;}
			else if (pm->bondInterType[nbond]==CS) { pfene = &fnCS;}
			else if (pm->bondInterType[nbond]==BB) { pfene = &fnBB;}
			else{
				printf("Wrong Fene\n");
			}
			fene += Fene(pfene,sqrt(rsq));
#endif
        }

    }
    double morse =0.0;
    for (monnext = 0; monnext < pd->numOfMonomers; monnext++){
        pmnext = &pd->monomer[monnext];
        if (pm!=pmnext){
            typem     = pm->type;
            typemnext = pmnext->type;

            dr.x =monpos->x-pmnext->pos.x;
            dr.y =monpos->y-pmnext->pos.y;
            dr.z =monpos->z-pmnext->pos.z;
#if USE_PBC==1
            PBCAll(dr);
#endif
            rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            
#ifdef LUT_MORSE
	 		if (typem == typemnext){
	 			if (typem == B || typemnext == C)     { interType = CC;}
	 			else if (typem == S || typemnext ==S) { interType = SS;}
	 		}else if (typem != typemnext){
	 			if (typem == S || typemnext == S)     { interType = CS;}
	 			else                                  { interType = CC;}
	 		}
	 		morse += LUTMorseValue(interType,rsq);
#else
			if (typem == typemnext){
				if (typem == B || typemnext == C)     { pmrs = &MrsCC;}
				else if (typem == S || typemnext ==S) { pmrs = &MrsSS;}
			}else if (typem != typemnext){
				if (typem == S || typemnext == S)     { pmrs = &MrsCS;}
				else                                  { pmrs = &MrsCC;}
			}
			if (rsq<SQR(pmrs->RCut)){
				morse += Morse(pmrs,sqrt(rsq))-Morse(pmrs,pmrs->RCut);
			}
#endif
        }
    }
    
    return fene+morse;
}
double EnergyMonomer (VecR pos,Monomer *pm,int mb){
// Calculates the energy of the mcur Monomer
// with monomers with id from mb to Nmon
// pos    : trial position of current monomer mcur
// typecur: type of current monomer Mon[mcur].type
    int    m,nbond,typecur,typenei;
    Monomer *pmnei;
    int indend;
    VecR   dr;
    double en;
    double r,r2;

    MorseParams* pmrs =NULL;
    FeneParams *pfene=NULL;

    en=0;
    typecur = pm->type;
    indend=pm->inDend;
    // FENE for bonded monomers
    for (nbond=0;nbond<pm->numOfBonds;nbond++){
        pmnei    = pm->bond[nbond];
        if (pmnei!=NULL){
            typenei =  pmnei->type;
            if (typecur == typenei){
                if      (typecur == B) { pfene = &fnBB;}
                else if (typecur == C) { pfene = &fnCC;}
            }
            if (typecur !=typenei){
                if (typecur == S || typenei == S) { pfene = &fnCS;}
                else                              { pfene = &fnCC;}
            }
            if ((pmnei->id>=mb)&&(pmnei!=NULL)){
                dr.x=pos.x-pmnei->pos.x;
                dr.y=pos.y-pmnei->pos.y;
                dr.z=pos.z-pmnei->pos.z;

                r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                r=sqrt(r2);
                en+=Fene(pfene,r);

            }
        }
    }
    // Calculate non bonded monomer interactions that are
    // not already calculated
    for (m=mb;m<dendrimer[indend].numOfMonomers;m++){
        pmnei = &dendrimer[indend].monomer[m];
        if (pm->id!=m){
            typenei = pmnei->type;
            if (typecur == typenei){
                if (typecur == B || typecur == C) {
                    pmrs = &MrsCC;

                }
                if (typecur == S){
                    pmrs = &MrsSS;

                }
            }
            if (typecur != typenei){
                if (typecur == S || typenei == S) {
                    pmrs = &MrsCS;

                }
                else{
                    pmrs = &MrsCC;
                }
            }
            dr.x=pos.x-pmnei->pos.x;
            dr.y=pos.y-pmnei->pos.y;
            dr.z=pos.z-pmnei->pos.z;
            r2=SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            if (r2<SQR(pmrs->RCut)){
                r=sqrt(r2);
                en += Morse(pmrs,r)-Morse(pmrs,pmrs->RCut);
            }
        }
    }
    return en;
}
double FeneEnergy(Monomer *pm,VecR *monpos,int mb){
    // mb: avoid calculating two times the same
    //     energy(call with monomer id +1 )
    int nbond;
    Dendrimer *pd;
    Monomer   *pb;
    VecR dr;
    double rsq,result;
#ifdef LUT_FENE
	int interType =-1;
#else
    FeneParams *pfene;
#endif
    result =0.0;
            
    pd = dendrimer + (pm->inDend);
    
    for (nbond = 0; nbond < pm->numOfBonds ; nbond++){
        pb = pm->bond[nbond];
        if ((pb != NULL)&&(pb->id>=mb)){
            dr.x = monpos->x-pd->monomer[pb->id].pos.x;
            dr.y = monpos->y-pd->monomer[pb->id].pos.y;
            dr.z = monpos->z-pd->monomer[pb->id].pos.z;
#if USE_PBC==1
            PBCAll(dr);
#endif
            rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
#ifdef LUT_FENE
            if      (pm->bondInterType[nbond]==CC) { interType = CC;}
            else if (pm->bondInterType[nbond]==CS) { interType = CS;}
            else if (pm->bondInterType[nbond]==BB) { interType = BB;}
			result += LUTFeneValue(interType,rsq);
#else
            if      (pm->bondInterType[nbond]==CC) { pfene = &fnCC;}
            else if (pm->bondInterType[nbond]==CS) { pfene = &fnCS;}
            else if (pm->bondInterType[nbond]==BB) { pfene = &fnBB;}
            else pfene = NULL;
            result += Fene(pfene,sqrt(rsq));            
#endif
        }
    }
    return result;
}
double EnergySystem(){
    // InterDend IntraDend
    double result =0.0;
#if USE_CELL_LIST==1
    double morse,fene;
    morse  = 0.0;
    fene   = 0.0;
    result = 0.0;
    struct strMonomer* pmon;
    for (int ndend=0; ndend< numOfDendrimers; ndend++){
        for (int mon=0;mon<dendrimer[ndend].numOfMonomers;mon++){
            pmon    = (dendrimer[ndend].monomer) + mon;
            fene   += FeneEnergy(pmon,&pmon->pos,0);
            morse  += CelldESingleMonMove(pmon,&pmon->pos);
        }
    }
    result     = fene/2.0 + morse/2.0;
//  
	MorseTotal = morse / 2.0;
	FeneTotal  = fene  / 2.0; 
#else
    result +=IntraDend(&dendrimer[numOfDendrimers-1]);
    for (int ndend=0; ndend< numOfDendrimers-1; ndend++){
        result += IntraDend(&dendrimer[ndend]);
        for (int ndendnext=ndend+1;ndendnext <numOfDendrimers;ndendnext++){
            result += InterDend((dendrimer+ndend),(dendrimer+ndendnext));
        }
    }
#endif
    return result;
}

double EnergySystem2(){
    double energy =0.0;
    Dendrimer *pd1,*pd2;

    for (int dendid1 =0; dendid1 < numOfDendrimers-1; dendid1++){
        pd1= &dendrimer[dendid1];
        energy += EnergyDend(pd1);
        for (int dendid2 =dendid1+1; dendid2 < numOfDendrimers; dendid2++){
                pd2 = &dendrimer[dendid2];
                energy +=PotDendDend(pd1,pd2);
        }
    }
    energy += EnergyDend(dendrimer+(numOfDendrimers-1));

    return energy;
}
double EnergySystem3(){
    double result = 0.0;
    for (int ndend=0; ndend< numOfDendrimers-1; ndend++){
        result += IntraDend(&dendrimer[ndend]);
        for (int ndendnext=ndend+1;ndendnext <numOfDendrimers;ndendnext++){
            result += InterDend((dendrimer+ndend),(dendrimer+ndendnext));
        }
    }
    result +=IntraDend(dendrimer+(numOfDendrimers-1));

    return result;
}

double EnergyDend (Dendrimer* pd){
    double energy = 0.0;
    for (int mon=0;mon < pd->numOfMonomers ; mon++){
        energy += PotMonIntraDend(&pd->monomer[mon],&pd->monomer[mon].pos,0);
    }
    energy/=2.0;
    return energy;
}

void ApplyBoundaryConditions(){
    int ndend, nmon;
    Dendrimer *pd;
    for (ndend=0;ndend<numOfDendrimers;ndend++){
        pd=dendrimer + ndend;
        for (nmon=0;nmon<pd->numOfMonomers;nmon++){
            PBCAll(pd->monomer[nmon].pos);
			//VWrapAll(pd->monomer[nmon].pos);
        }
    }
}
void CalcInterForce(Dendrimer* d1, Dendrimer* d2,VecR *forceCM1,VecR *forceCM2){
//  Calculates the total force acting from  molecule d2
//  to molecule d1 . The force can be found by integrating the
//  Morse potential energy (MorceForce () )
//  The force is then projected to the line connecting the centre of masses
	Monomer *pm1,*pm2;
	int mon1,mon2;
	int typem1,typem2;
	double force;
	VecR dr;
	VecR vCM01,vCM10;
	double norma01;
	double dummyProj;

	MorseParams* pmrs=NULL;

	forceCM1->x+=0.0;
	forceCM1->y+=0.0;
	forceCM1->z+=0.0;

	forceCM2->x-=0.0;
	forceCM2->y-=0.0;
	forceCM2->z-=0.0;

	for (mon1=0;mon1<d1->numOfMonomers;mon1++){
		pm1 = (d1->monomer)+mon1;
		typem1 = pm1->type;
		for (mon2=0;mon2<d2->numOfMonomers;mon2++){
			pm2 = (d2->monomer)+mon2;
			typem2 = pm2->type;
			if (typem1 == typem2){
				if (typem1 == B || typem2 == C) pmrs = &MrsCC;
				else if (typem1 == S)           pmrs = &MrsSS;
			}else if (typem1 != typem2){
				if (typem1 == S || typem2 == S) pmrs = &MrsCS;
				else                            pmrs = &MrsCC;
			}
			dr.x = (pm1->pos).x - (pm2->pos).x; 
			dr.y = (pm1->pos).y - (pm2->pos).y; 
			dr.z = (pm1->pos).z - (pm2->pos).z; 
			double dumnorma = VLen(dr);
			force=MorseForce(pmrs,dr); // force on mon1 from mon2

			forceCM1->x+=force*dr.x/dumnorma;
			forceCM1->y+=force*dr.y/dumnorma;
			forceCM1->z+=force*dr.z/dumnorma;

			forceCM2->x-=force*dr.x/dumnorma;
			forceCM2->y-=force*dr.y/dumnorma;
			forceCM2->z-=force*dr.z/dumnorma;
		}
	}
	DendCalcCMass(d1);
	DendCalcCMass(d2);
	
	vCM01.x = d1->CMass.x - d2->CMass.x;
	vCM01.y = d1->CMass.y - d2->CMass.y;
	vCM01.z = d1->CMass.z - d2->CMass.z;
	norma01=VLen(vCM01);
	vCM01.x = vCM01.x/norma01;
	vCM01.y = vCM01.y/norma01;
	vCM01.z = vCM01.z/norma01; 	
	vCM10.x = -vCM01.x;
	vCM10.y = -vCM01.y;
	vCM10.z = -vCM01.z;
	
	dummyProj = VDot((*forceCM1),vCM01)/VDot(vCM01,vCM01);
	forceCM1->x = dummyProj*forceCM1->x;
	forceCM1->y = dummyProj*forceCM1->y;
	forceCM1->z = dummyProj*forceCM1->z;

	dummyProj = VDot((*forceCM2),vCM10)/VDot(vCM10,vCM10);
	forceCM1->x = dummyProj*forceCM1->x;
	forceCM1->y = dummyProj*forceCM1->y;
	forceCM1->z = dummyProj*forceCM1->z;

}

double  MorseForce(MorseParams* pmrs,VecR dr){
// calculates the force vector between two monomers 
// INPUT:  distance vector dr
// OUTPUT: force vector f 
	double eps,a,b,d;
	double force,rr;
	eps    = pmrs->eps;
	a      = pmrs->a;
	b      = 1.0;
	d      = pmrs->d;
	rr     = sqrt(SQR(dr.x)+SQR(dr.y)+SQR(dr.z));
	if (rr<pmrs->RCut) force  = (2.0*eps*a/b) * exp(-a*(rr-d) )*( exp(-a*(rr-d)) -1 );
	else force=0.0; 

	return force;
}
// usefull function for testing purposes

double MorseEnergy(Monomer *pm, VecR *monpos){
    Dendrimer *pd;
    int typem;
    int monnext, typemnext;
    Monomer *pmnext;
    int mb=0;
    VecR dr;
    double rsq=0.0,result=0.0;
    pd = &dendrimer[pm->inDend];

    MorseParams *pmrs=NULL;

    for (monnext = mb; monnext<pd->numOfMonomers; monnext++){
        if (pm->id!=monnext){
            pmnext = &pd->monomer[monnext];
            typem     = pm->type;
            typemnext = pmnext->type;

            dr.x =monpos->x-pmnext->pos.x;
            dr.y =monpos->y-pmnext->pos.y;
            dr.z =monpos->z-pmnext->pos.z;
#if USE_PBC==1
            PBCAll(dr);
#endif
            rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
            
			if (typem == typemnext){
				if (typem == B || typemnext == C) pmrs = &MrsCC;
				else if (typem == S)              pmrs = &MrsSS;
			}else if (typem != typemnext){
				if (typem == S || typemnext == S) pmrs = &MrsCS;
				else                              pmrs = &MrsCC;
			}
			if (rsq<SQR(pmrs->RCut)){
				result += Morse(pmrs,sqrt(rsq))-Morse(pmrs,pmrs->RCut);
			}
        }
    }
    return result;
}

void EnergyCheck(){
    // check energy
    double energysystem1 ;
    double energysystem2 ;
    double energysystem3 ;
    energysystem1 = EnergySystem();
    energysystem2 = EnergySystem2();
    energysystem3 = EnergySystem3();
    printf("Total Energy System->               "
            " \n1  :%30.14f \n2  :%30.14f \n3  :%30.14f \nRE :%30.14f\n",
            energysystem1,energysystem2,energysystem3,RunningEnergy);
	fflush(stdout);

#if USE_CELL_LIST==10   
    Monomer *pm;
    double cellenergy, cellenergy2, totalenergy1, totalenergy2 ;
    double bondedenergy=0.0;
    cellenergy = 0.0; cellenergy2=0.0;
    for (int i=0;i<numOfDendrimers;i++){
        for (int j=0;j<dendrimer[i].numOfMonomers;j++){
            pm = &dendrimer[i].monomer[j];
            bondedenergy+=FeneEnergy(pm,&pm->pos,0);
            cellenergy2+=CelldESingleMonMove(pm,&pm->pos);
        }
    }
    
    cellenergy   = CellInterDend();
    totalenergy1 = bondedenergy/2.0 + cellenergy;
    totalenergy2 = bondedenergy/2.0 + cellenergy2/2.0;

    if (IsEqual(totalenergy1,totalenergy2,1E-7)==FALSE){
        printf("Error: EnergyCheck!! Energy different, %lf %lf\n",totalenergy1,totalenergy2);
        printf("Cell Energy Values:%lf,%lf\n",totalenergy1, totalenergy2);
        printf("exit program\n");
    }
#endif

}

int IsEqual(double a, double b, double tol){
    if (fabs(a-b)<tol){
        return TRUE;
    }else{
        return FALSE;
    }
}
