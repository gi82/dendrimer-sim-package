/* 
 * File:   iniconfig.c
 * Author: georgiou
 *
 * Created on May 19, 2011, 2:18 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include "system.h"
/*
 * 
 */
int main(int argc, char** argv) {
    int numOfInParams = 2;
    if (argc<numOfInParams){
        printf("two few arguments while calling %-s:\n",argv[0]);
        printf("number of Params:%-2d\n",argc-1);
        exit(EXIT_FAILURE);
    }
    else if (argc==numOfInParams){
        ReadInputData(argv[1]);
    }
    
    return (EXIT_SUCCESS);
}