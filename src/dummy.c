#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define VWrap(v, t) \
if (v.t >= 0.5 * Box.t) v.t -= Box.t; \
else if (v.t < -0.5 * Box.t) v.t += Box.t
#define VWrapAll(v) \
{VWrap (v, x); \
VWrap (v, y);\
VWrap (v, z);}



typedef struct {
    double x,y,z;
}  VecR ;

typedef struct{
    int x,y,z;
} VecI;
#define PBC1_New(v,v0,t) ( (v).t = (v).t - Box.t *lround(  ( (v).t - (v0).t ) / Box.t ) )

#define PBC1_halb(v,t) \
    if      ( (v).t >  0.5*Box.t ) { (v).t -= (int)( 2*(v).t /Box.t )*Box.t;} \
    else if ( (v).t < -0.5*Box.t ) { (v).t += (int)( 2*(v).t /Box.t )*Box.t;}
// PBC (0..Box)
#define PBC1_BOX(v,t) \
    if      ( (v).t > Box.t) {(v).t -= (int)( (v).t /Box.t )*Box.t; } \
    else if ( (v).t < 0.0  ) {(v).t += (1-(int)( (v).t /Box.t ))*Box.t;}
// PBC (-Box/2.0 ... Box/2.0)
#define PBC1(v,t) ( (v).t = (v).t - lround(((v).t)/Box.t)*(Box.t) )
#define PBCAll(v)   \
{   PBC1(v,x);      \
    PBC1(v,y);      \
    PBC1(v,z); }

VecI Box;

static void DendPBCunwrap(VecR* pos, const VecR *refpos);
int main(){
    
    Box.x =20;
    Box.y=20;
    Box.z=20;

    VecR test;
    VecR ref={3.0,3.0,3.0};

    fprintf(stdout,"Box size: %d, %d, %d\n",Box.x, Box.y, Box.z);
    for (;;){
        fprintf(stdout,"Give vector x y z\n");
        scanf("%lf %lf %lf",&test.x,&test.y,&test.z);
        fprintf(stdout,"test vector : %lf %lf %lf\n",test.x, test.y, test.z);
        fprintf(stdout,"Wrapping\n");
        //PBCAll(test);
        PBCAll(test);
        fprintf(stdout,"new test vector : %lf %lf %lf\n",test.x, test.y, test.z);
    }

    return 0;
}
static void DendPBCunwrap(VecR* pos, const VecR *refpos){
    //pos->x -= lround(  ( pos->x - refpos->x ) / Box.x )*(Box.x);
    //pos->y -= lround(  ( pos->y - refpos->y ) / Box.y )*(Box.y);
    //pos->z -= lround(  ( pos->z - refpos->z ) / Box.z )*(Box.z);

    pos->x -= lround( ( pos->x - refpos->x ) / Box.x )*(Box.x);
    pos->y -= lround( ( pos->y - refpos->y ) / Box.y )*(Box.y);
    pos->z -= lround( ( pos->z - refpos->z ) / Box.z )*(Box.z);
    printf("lround values: %lf, %ld, %lf, %ld, %lf, %ld\n"
            ,pos->x - refpos->x,lround(pos->x - refpos->x)
            ,pos->y - refpos->y,lround(pos->y - refpos->y)
            ,pos->z - refpos->z,lround(pos->z - refpos->z)
            );
}
/*
#define PBC1(v,t) ( (v).t = (v).t - lround(((v).t)/Box.t)*(Box.t) )
double  P_Cd_x (double x){
	double x1;
	x1=x;
	if (x>XBOX) x1=x-(int)(x/XBOX)*XBOX;

	if (x<0)	x1=x+(1-(int)(x/XBOX))*XBOX;

	return (x1);
}
double  P_Cd_y (double y){
	 double y1;
	y1=y;
	if (y>YBOX) y1=y-(int)(y/YBOX)*YBOX;

	if (y<0)	y1=y+(1-(int)(y/YBOX))*YBOX;

	return (y1);
}
double  P_Cd_z (double z){
	 double z1;
	z1=z;
	if (z>ZBOX) z1=z-(int)(z/ZBOX)*ZBOX;

	if (z<0)	z1=z+(1-(int)(z/ZBOX))*ZBOX;

	return (z1);
}
*/
