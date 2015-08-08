
#ifndef FBCELLFINDER_SUBS_H
#define FBCELLFINDER_SUBS_H

#include <stdio.h>




void cutfname(char *infile, char *galID, char *ion, char *lostag, char *outfile);

void gethdr( int klos, double a, double *xg, double *yg, double *zg, double *vxg, double *vyg, double *vzg, double *b1, double *b2, double *x0, double *y0, double *z0, double *l, double *m, double *n, double ap[][3], double *zbox, double *vgal, double *zgal, char *infile);

void readcells( int *cellnum, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *Lcell, double *ndencell, double *fion, double *temp, double *zmfrac, char *infile);

#endif

