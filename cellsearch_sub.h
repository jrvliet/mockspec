#ifndef CELLSEARCH_SUBS_H
#define CELLSEARCH_SUBS_H

#include <stdio.h>

void cellsearch(int nlos, char *GZfile, double *x0, double *y0, double *z0, double *lb, double *mb, double *nb, double *xen, double *yen, double *zen, double *xex, double *yex, double *zex, double *lowervel, double *uppervel, FILE *fpout, FILE *fplog, double *a11, double *a12, double *a13, double *a21, double *a22, double *a23, double *a31, double *a32, double *a33);

double get_dmin(double cellsize, double x, double y, double z, double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z);


void confirm_cell(double cellsize, double x, double y, double z, double dmin, double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z, double Dx, double Dy, double Dz, double r_max, double x0, double y0, double z0, int *in_cell, int *onWall, double *rx_xpt, double *ry_xpt, double *rz_xpt);


double get_radvel(double x, double y, double z, double vx, double vy, double vz);



#endif
