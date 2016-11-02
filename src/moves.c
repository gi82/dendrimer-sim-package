#include "system.h"
#include "ran250.h"
// move centre of mass of dendrimer dendid
// move all monomers by vecr.x,y,z
void DendPlaceCMass(Dendrimer *pd,VecR vecr){
// USES npos
// puts centre of mass of dendrimer to
// specific position
    VecR tvecr;
    tvecr.x = -pd->CMass.x;
    tvecr.y = -pd->CMass.y;
    tvecr.z = -pd->CMass.z;

    DendMoveCMass(pd,tvecr);
    DendMoveCMass(pd,vecr);

}
void DendMoveCMass(Dendrimer *pd, VecR vecr) {
    int first, last;
    first = 0; //dendrimer[dendid].first;
    last = pd->numOfMonomers; //dendrimer[dendid] .last;
    for (int mon = first; mon < last; mon++) {
        pd->monomer[mon].npos.x += vecr.x;
        pd->monomer[mon].npos.y += vecr.y;
        pd->monomer[mon].npos.z += vecr.z;
    }
    pd->CMass.x += vecr.x;
    pd->CMass.y += vecr.y;
    pd->CMass.z += vecr.z;
}
// Move monomer monid and then restore CMass position
// by sybstracting vecr/numOfMonomer from the rest of
// the monomers in the same dendrimer.

void MonSingleMove(Monomer *pm, VecR* pvecr, int FixedCMass) {
    int first, last;

    pm->bupos.x = pm->pos.x;
    pm->bupos.x = pm->pos.x;
    pm->bupos.x = pm->pos.x;

    pm->pos.x += (*pvecr).x;
    pm->pos.y += (*pvecr).y;
    pm->pos.z += (*pvecr).z;
    if (FixedCMass == TRUE) {
        Dendrimer *pd = &dendrimer[pm->inDend];
        first = 0;
        last = pd->numOfMonomers;
        for (int mon = first; mon < last; mon++) {
            if (mon != pm->id) {
                pm->pos.x -= pvecr->x / (double) (pd->numOfMonomers - 1);
                pm->pos.y -= pvecr->y / (double) (pd->numOfMonomers - 1);
                pm->pos.z -= pvecr->z / (double) (pd->numOfMonomers - 1);
            }
        }
    }
}
void MonTrialMove_All(VecR *pnewPBCpos, const VecR *poldPBCpos,
                      VecR* pnewNpos,   const VecR *poldNpos) {

    double radius;
    VecR vecr;
    radius = MaxStep * dr250();

    RandomVector(&vecr, radius);
//  update pbc position
    pnewPBCpos->x = poldPBCpos->x + vecr.x;
    pnewPBCpos->y = poldPBCpos->y + vecr.y;
    pnewPBCpos->z = poldPBCpos->z + vecr.z;
//  update non pbc position
    pnewNpos->x   = poldNpos->x   + vecr.x;
    pnewNpos->y   = poldNpos->y   + vecr.y;
    pnewNpos->z   = poldNpos->z   + vecr.z;

}
void MonTrialMove(VecR *pnewPos, const VecR* poldPos) {
    double radius;
    VecR vecr;
    radius = MaxStep * dr250();

    RandomVector(&vecr, radius);

    pnewPos->x = poldPos->x + vecr.x;
    pnewPos->y = poldPos->y + vecr.y;
    pnewPos->z = poldPos->z + vecr.z;
}
void MonTrialMoveFixedCMass(Monomer *pm){
    
    Dendrimer *pd;
    Monomer   *pmnext;
    double radius;
    VecR vecr;
    
    pd = &dendrimer[pm->inDend];

    pm->bupos.x = pm->pos.x;
    pm->bupos.y = pm->pos.y;
    pm->bupos.z = pm->pos.z;
    
    radius = MaxStep * dr250();
    RandomVector(&vecr, radius);
    
    pm->pos.x += vecr.x;
    pm->pos.y += vecr.y;
    pm->pos.z += vecr.z;
    for (int mon = 0; mon < pd->numOfMonomers; mon++) {
        pmnext = &pd->monomer[mon];
        if (pm->id!=pmnext->id){
            pmnext->bupos.x = pmnext->pos.x;
            pmnext->bupos.y = pmnext->pos.y;
            pmnext->bupos.z = pmnext->pos.z;

            pmnext->pos.x -= vecr.x / (double) (pd->numOfMonomers-1);
            pmnext->pos.y -= vecr.y / (double) (pd->numOfMonomers-1);
            pmnext->pos.z -= vecr.z / (double) (pd->numOfMonomers-1);
        }
    }
}

void RandomVector(VecR* pvecr, double radius) {
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
    pvecr->z = radius * (1.0 - 2.0 * ransq);

    ranh = 2.0 * sqrt(1 - ransq);
    radius = radius*ranh;

    pvecr->x = radius*ran1;
    pvecr->y = radius*ran2;

}

void WriteVecR(VecR vecr, FILE *pOutFile) {
    fprintf(pOutFile, "%12.6f  %12.6f  %12.6f\n", (vecr).x, (vecr).y, (vecr).z);
}
void WriteVecI(VecI veci, FILE *pOutFile) {
    fprintf(pOutFile, "%d  %d  %d\n", veci.x, veci.y, veci.z);
}

/*void MonTrialMoveTwoFixedCMass(Monomer *pm1,Monomer *pm2){
    double radius;
    VecR   vecr;
    // moves monomers pm1 & pm2 to opposite directions
    // so as the Centre of mass of the dendrimers remains
    // fixed
    // pm1 & pm2 have to be on the SAME DENDRIMER
    radius = MaxStep * dr250();
    RandomVector(&vecr, radius);

    pm1->bupos.x = pm1->pos.x;
    pm1->bupos.y = pm1->pos.y;
    pm1->bupos.z = pm1->pos.z;

    pm2->bupos.x = pm2->pos.x;
    pm2->bupos.y = pm2->pos.y;
    pm2->bupos.z = pm2->pos.z;

    pm1->pos.x   += vecr.x;
    pm1->pos.y   += vecr.y;
    pm1->pos.z   += vecr.z;

    pm2->pos.x   -= vecr.x;
    pm2->pos.y   -= vecr.y;
    pm2->pos.z   -= vecr.z;

}*/


