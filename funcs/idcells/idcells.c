
/* 
    Reads in the list of celsl found by cellfinder.c
    Finds the cells in the ion boxes
    Passes it to los7
*/

#include "idcells-subs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int idcells(char *ion, char *galID, char *expn, int numcores, int losnum, 
            double impact, double phi){

    int i, losnum;
    double losB, losPhi;
    char boxfile[100];
    // Read in the LOS properties, stored in lines.info file
    
    // Read in the gas box
    buildBoxName(galID, expn, ion, boxfile);
    
    // Get the number of cells in the box
    numCells = boxSize(boxfile);

    // Allocate memory for the box file
    double *size = (double *)calloc(numCells, sizeof(double)
    double *x = (double *)calloc(numCells, sizeof(double)
    double *y = (double *)calloc(numCells, sizeof(double)
    double *z = (double *)calloc(numCells, sizeof(double)
    double *vx = (double *)calloc(numCells, sizeof(double)
    double *vy = (double *)calloc(numCells, sizeof(double)
    double *vz = (double *)calloc(numCells, sizeof(double)
    double *dense = (double *)calloc(numCells, sizeof(double)
    double *temp = (double *)calloc(numCells, sizeof(double)
    double *snII = (double *)calloc(numCells, sizeof(double)
    double *snIa = (double *)calloc(numCells, sizeof(double)
    
    
    FILE *fp = fopen(boxfile, "r");
    while
    







    

    return 0;
}

