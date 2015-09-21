
#ifndef CELLPATHSUB7_H
#define CELLPATHSUB7_H

#include <stdio.h>


void getD(double l, double m, double n, double x0, double y0, double z0, 
          double X, double Y, double Z, double lengthcell, double *cellpath, 
          int *error, int *errtype, int losnum);

void printArr( double *x, int size);



#endif
