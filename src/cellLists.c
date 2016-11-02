#include "system.h"

// call CellListSetup at initialization step after Box is defined
static double myrint(double x) {
  return ( (x<0.0)? -floor(-x+0.5):floor(x+0.5) );
}
static  int CellPeriodic1(int x,int NCells){
    int xx;
    xx=x;
    if (x>=NCells) xx = x - myrint(x/NCells)*NCells;
    if (x<0)	   xx = x + (1-myrint((x+1)/NCells))*NCells;

    return (xx);
}
static int CellPeriodic(int ix, int iy, int iz, VecI cellni){
 // Old Version
/*
    if (ix<0) ix += cellni.x;
    if (iy<0) iy += cellni.y;
    if (iz<0) iz += cellni.z;

    if (ix>=cellni.x) ix -= cellni.x;
    if (iy>=cellni.y) iy -= cellni.y;
    if (iz>=cellni.z) iz -= cellni.z;
*/
    ix = CellPeriodic1(ix,cellni.x);
    iy = CellPeriodic1(iy,cellni.y);
    iz = CellPeriodic1(iz,cellni.z);
    return (CELLINDX(ix,iy,iz,cellni));
}
void CellListSetup(VecI vBox,double rangemax){
    // using floor for bigger cells
    if (rangemax<=0.0){
        fprintf(stderr,"cellLists.c:: CellListSetupfailed, rc:%lf<0.0\n",rangemax);
        exit(4);
    }
    cellNI.y   = (int) floor((double)vBox.y/rangemax);
    cellNI.x   = (int) floor((double)vBox.x/rangemax);
    cellNI.z   = (int) floor((double)vBox.z/rangemax);
    
    cellSizeR.x   = (double)vBox.x / (double)cellNI.x;
    cellSizeR.y   = (double)vBox.y / (double)cellNI.y;
    cellSizeR.z   = (double)vBox.z / (double)cellNI.z;

    if ((cellSizeR.x<rangemax)||(cellSizeR.y<rangemax)||(cellSizeR.z<rangemax)){
        fprintf(stderr,"ERROR: length of subbox %g %g %gsmaller than "
                "interaction range %g\n",cellSizeR.x,cellSizeR.y,cellSizeR.z,rangemax);
        exit(4);
    }

    cellInvSizeR.x = 1.0/cellSizeR.x;
    cellInvSizeR.y = 1.0/cellSizeR.y;
    cellInvSizeR.z = 1.0/cellSizeR.z;

    // NumberofTotal cells in the Box
    //cellNYZ         = cellNI.y * cellNI.z;
    cellTotalNumber = cellNI.x * cellNI.y * cellNI.z;
}
// call CellAlloc() at Allocate step
void CellAlloc(){
    cells = (Cell*) malloc(cellTotalNumber*sizeof(Cell));
}
void CellDeAlloc(){
	if (cells!=NULL) free(cells);
	else printf("Cell List deallocate failed... Cell list is empty!!\n");
}
// Called once after all dendrinmers' positions
// have been defined
void CellNewBuild(){
    int  cx,cy,cz,cindex;
    int  neix, neiy, neiz, nei;
    int  dendid;
    int  mon;
    Cell *pcell;
    VecI vci;

    Dendrimer *pd;
    Monomer   *pm;
#if DEBUG_CELL_LIST==1
    printf("  %s:\n","Building New CellList!");
    printf("Cell List Size xyz:%3d %3d %3d\n",cellNI.x,cellNI.y,cellNI.z);
#endif
    // first Apply boundary conditions
#if USE_PBC==1
    ApplyBoundaryConditions();
#endif
    // Initialize all cells
    for (cx=0;cx<cellNI.x;cx++)
        for (cy=0;cy<cellNI.y;cy++)
            for (cz=0;cz<cellNI.z;cz++){
                cindex = CELLINDX(cx,cy,cz,cellNI);
#if DEBUG_CELL_LIST==1
                printf("cindex:%d,(%d,%d,%d)\n",cindex,cx,cy,cz);
#endif
                pcell = cells + cindex;
                pcell->numOfElements = 0;
                pcell->numOfNeigh    = 0;
                pcell->first = NULL;
                pcell->index = -1;
                // first neighbours around cell including itself
                nei = 0;
                for (neix=-1;neix<=1;neix++)
                    for (neiy=-1;neiy<=1;neiy++)
                        for (neiz=-1;neiz<=1;neiz++){
                            (pcell->neighbours)[nei]=CellPeriodic(cx+neix,cy+neiy,cz+neiz,cellNI);
                            nei++;
                }
                pcell->numOfNeigh = nei;
                // second neighbours around cells 5*5*5
                for (neix=-2;neix<=2;neix++)
                    for (neiy=-2;neiy<=2;neiy++)
                        for (neiz=-2;neiz<=2;neiz++){
                            if ( (abs(neix)==2) || (abs(neiy)==2) || (abs(neiz)==2) ){
                                (pcell->neighbours)[nei]=CellPeriodic(cx+neix,cy+neiy,cz+neiz,cellNI);
                                nei++;
                            }
               }
#if CELL_SEC_NEI==1
                pcell->numOfNeigh = nei;
#endif
    }

    
    for (int ndend=0;ndend<numOfDendrimers;ndend++){
        pd = &dendrimer[ndend];
        for (int nmon=0;nmon<pd->numOfMonomers;nmon++){
            pd->monomer[nmon].previous = NULL;
            pd->monomer[nmon].next = NULL;
        }
    }
    for (cindex=0;cindex<cellTotalNumber;cindex++){
        (cells+cindex)->index = cindex;
        (cells+cindex)->first = NULL;
        (cells+cindex)->numOfElements = 0;
    }

    printf("putting monomers in cells\n");
    for (dendid=0;dendid<numOfDendrimers;dendid++){
        pd = &dendrimer[dendid];
        for (mon=0;mon<pd->numOfMonomers;mon++){
            pm = (pd->monomer)+mon;
            
            vci.x = IFLOAT(pm->pos.x,cellSizeR.x);
            vci.y = IFLOAT(pm->pos.y,cellSizeR.y);
            vci.z = IFLOAT(pm->pos.z,cellSizeR.z);

            cindex = CellPeriodic(vci.x,vci.y,vci.z,cellNI);
#if DEBUG_CELL_LIST==1
            printf("cindex:%d,(%lf,%lf,%lf)\n",cindex,pm->pos.x,pm->pos.y,pm->pos.z);
#endif
            CellAddElement(cindex,pm);
        }
    }
    int count;
    for (count=0,cindex=0;cindex<cellTotalNumber;cindex++){
        pcell = cells + cindex;
        count+=pcell->numOfElements;
    }
    if (count!=numOfTotalMonomers) {
        printf("!!! not all monomers are in cells\n");
    }
    else{
        printf("All monomers in cells\n");
        WriteVecI(cellNI,stdout);
        WriteVecR(cellSizeR,stdout);
    }
}
// Put monomer pmon to Cell(index)
void CellUpdateMon(struct strMonomer *pmon){
    int cindex_old,cindex_new;
    VecI vci;
    cindex_old = pmon->cellIndex;
    
    vci.y = IFLOAT(pmon->pos.y,cellSizeR.y);
    vci.x = IFLOAT(pmon->pos.x,cellSizeR.x);
    vci.z = IFLOAT(pmon->pos.z,cellSizeR.z);
    cindex_new = CellPeriodic(vci.x,vci.y,vci.z,cellNI);

    if (cindex_new!=cindex_old){

        CellDelElement(cindex_old,pmon);
        CellAddElement(cindex_new,pmon);
    }
}
void CellUpdateAll(struct strDendrimer *pd){
    int mon;
    struct strMonomer* pm;
    for (mon=0;mon<pd->numOfMonomers;mon++){
        pm=(pd->monomer)+mon;
        CellUpdateMon(pm);
    }
}
int CellUpdateList(){
    int dendid;
    for (dendid=0; dendid<numOfDendrimers; dendid++){
        CellUpdateAll( (dendrimer + dendid) );
    }
    return 0;
}

void CellAddElement(int cellindex, struct strMonomer *pmon){
    struct strMonomer *pprevious;
    Cell *pcell=cells+cellindex;// &cells[cellindex]
// if cell is not empty
    if (pcell->numOfElements != 0){
        pprevious                = (pcell->first)->previous;
        pprevious->next          = pmon;
        (pcell->first)->previous = pmon;
        pmon->previous           = pprevious;
        pmon->next               = pcell->first;
        pmon->cellIndex          = cellindex;
        pcell->numOfElements++;

    }else{
        pmon->cellIndex          = cellindex;
        pcell->first             = pmon;
        pmon->previous           = pmon;
        pmon->next               = pmon;
        pcell->numOfElements++;
    }

}

void CellDelElement(int cellindex, struct strMonomer *pmon){
   
    struct strMonomer *pprevious,*pnext;
    Cell *pcell = cells + cellindex;

    pprevious = pmon->previous;
    pnext     = pmon->next;

    pprevious->next = pnext;
    pnext->previous = pprevious;
    pcell->first    = pprevious;
    pcell->numOfElements--;
    pmon->cellIndex = -100;

}

void CellTest(){
    int  cindex;
	//VecI cellI;
    FILE *pOutFile;
    pOutFile = stdout;

    fprintf(pOutFile,"Number of cells in x,y,z dir:");
    WriteVecI(cellNI,pOutFile);
    fprintf(pOutFile,"Size in x,y,z dir:           ");
    WriteVecR(cellSizeR,pOutFile);
    fprintf(pOutFile,"inverse Size in x,y,z dir   :");
    WriteVecR(cellInvSizeR,pOutFile);
    fprintf(pOutFile,"CELL LIST CONTENT\n");
    for (cindex=0;cindex<cellTotalNumber;cindex++){
        CellOutElements(cindex);
        //CELLCOORDS(cellI,cindex,cellNI);
    }

}

double CellInterDend(){
    // 
    double   result=0.0;
    int      cindex, nei;
    Cell    *pcell1,*pcell2;
    Monomer *pmon1 ,*pmon2;
    int      typem1, typem2;
    VecR     dr;
    double   rsq;
    int      mon1, mon2;

    MorseParams *pmrs=NULL;

    // Loop over all cells
    for (cindex= 0;cindex<cellTotalNumber;cindex++){
        pcell1 = cells + cindex;
        pmon1  = pcell1->first;
        for (mon1=0;mon1<pcell1->numOfElements;mon1++){
            for (nei=0;nei<pcell1->numOfNeigh;nei++){
                pcell2 = cells + (pcell1->neighbours)[nei];
                pmon2  = pcell2->first;
                // loop over elements of cells and neighbouring cells
                    for (mon2=0;mon2<pcell2->numOfElements;mon2++){
                        if ((pmon1!=NULL)&&(pmon2!=NULL)){
                            if (pmon1!= pmon2){
                                typem1=pmon1->type;
                                typem2=pmon2->type;
                                if (typem1 == typem2){
                                    if (typem1 == B || typem2 == C) pmrs = &MrsCC;
                                    else if (typem1 == S)           pmrs = &MrsSS;
                                }else if (typem1 != typem2){
                                    if (typem1 == S || typem2 == S) pmrs = &MrsCS;
                                    else                            pmrs = &MrsCC;
                                }
                                dr.x = (pmon1->pos).x - (pmon2->pos).x;
                                dr.y = (pmon1->pos).y - (pmon2->pos).y;
                                dr.z = (pmon1->pos).z - (pmon2->pos).z;
                                PBCAll(dr);
                                rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                                result += Morse(pmrs,sqrt(rsq))-Morse(pmrs,pmrs->RCut);
                            }
                        }
                        pmon2 = pmon2->next;
                    }

                }
                pmon1 = pmon1->next;
        }
    }
    return result/2.0;
    
}


void CellOutElements(int cellindex){
    Cell *pcell;
    Monomer *pmon;
    VecI cellI;
    int nei, cindex;

    pcell = cells + cellindex;
    pmon = pcell->first;
    CELLCOORDS(cellI,cellindex,cellNI);
    printf("cellnr:%d",pcell->index);
    printf("(%3d %3d %3d)\n",cellI.x,cellI.y,cellI.z);
    for (int elem=0;elem<pcell->numOfElements;elem++){
        printf(" inDend:%3d, monid:%3d, cellindx:%3d\n",pmon->inDend,pmon->id,pmon->cellIndex);
        pmon = pmon->next;
    }
    printf("neighbourhring cells:");
    for (nei=0;nei<pcell->numOfNeigh;nei++){
        cindex = (pcell->neighbours)[nei];
        CELLCOORDS(cellI,cindex,cellNI);
        //printf("%3d(%3d %3d %3d),",cindex,cellI.x,cellI.y,cellI.z);
    }
    printf("\n");

}


/*double CellInterDend1(){
    double   result=0.0;
    int      cindex, nei;
    Cell    *pcell1,*pcell2;
    Monomer *pmon1 ,*pmon2;
    int      typem1, typem2;
    VecR     dr;
    double   rsq;
    int      mon1, mon2;

    MorseParams *pmrs;
    // Loop over all cells
    for (cindex= 0;cindex<cellTotalNumber;cindex++){
        pcell1 = cells + cindex;
        pmon1  = pcell1->first;
        for (nei=0;nei<pcell1->numOfNeigh;nei++){
            pcell2 = cells + (pcell1->neighbours)[nei];
            pmon2  = pcell2->first;
            // loop over elements of cells and neighbouring cells
            if ((pmon1!=NULL)&&(pmon2!=NULL)){
                for (mon1=0;mon1<pcell1->numOfElements;mon1++){
                    for (mon2=0;mon2<pcell2->numOfElements;mon2++){
                        if (pmon1!= pmon2){
                            typem1=pmon1->type;
                            typem2=pmon2->type;
                            if (typem1 == typem2){
                                if (typem1 == B || typem2 == C) pmrs = &MrsCC;
                                else if (typem1 == S)           pmrs = &MrsSS;
                            }else if (typem1 != typem2){
                                if (typem1 == S || typem2 == S) pmrs = &MrsCS;
                                else                            pmrs = &MrsCC;
                            }
                            dr.x = (pmon1->pos).x - (pmon2->pos).x;
                            dr.y = (pmon1->pos).y - (pmon2->pos).y;
                            dr.z = (pmon1->pos).z - (pmon2->pos).z;
                            PBCAll(dr);
                            rsq  = SQR(dr.x)+SQR(dr.y)+SQR(dr.z);
                            result += Morse(pmrs,sqrt(rsq))-Morse(pmrs,pmrs->RCut);
                        }
                        pmon2 = pmon2->next;
                    }
                    pmon1 = pmon1->next;
                }

            }
        }
    }
    return result/2.0;
}*/
