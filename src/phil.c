
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
        	if (MaxStep/curstep>Box.x/4.0) MaxStep = Box.x/4.0;
		}

    }

}
