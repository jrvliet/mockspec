#ifndef GENLOS_FUNCS_H
#define GENLOS_FUNCS_H

#include <stdio.h>

int inBox(double x, double y, double z, double boxsize);
void findEnds(double px, double py, double pz, double dx, double dy, double dz, double boxsize, double *xen, double *yen, double *zen, double *xex, double *yex, double *zex);
void read_summary(char *galID, double *aexpn, char *summaryLoc, double *mvir, double *rvir, double *a11, double *a12, double *a13, double *a21, double *a22, double *a23, double *a31, double *a32, double *a33, double *vpec_x, double *vpec_y, double *vpec_z);
double toRadians( double degrees );




#endif
