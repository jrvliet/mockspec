
#ifndef FBCELLFINDER_SUBS_H
#define FBCELLFINDER_SUBS_H

#include <stdio.h>




void cutfname(char *infile, char *galID, char *ion, char *lostag, char *outfile);

void gethdr( int klos, double *a, double *xg, double *yg, double *zg, double *vxg, double *vyg, double *vzg, double *b1, double *b2, double *x0, double *y0, double *z0, double *l, double *m, double *n, double ap[][3], double *zbox, double *vgal, double *zgal, char *infile);

int readcells( int *cellnum, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *Lcell, double *ndencell, double *fion, double *temp, double *zmfrac, char *infile);

double getamu(char *tranilist, char *ionlabel);

void mkfname(char *infile, char *galID, char *ion, char *lostag, char *linesfile);

void wrtlines(double zgal, double *zline, double *Nline, double *bline, int *cellnum, char *linesfile, int ndata);

void wrtlosdata( double Slos, double Rgal, double zline, double vlos, double vabs, double dlos, double ndencell, double fion, double zmfrac, double Nline, double temp, double bline, double Vgalt, double vrp, double V_theta, double V_phi, double vzp, double xp, double yp, double zp, double rp, double theta, double phi, int cellnum, char *unitlosfile);

void deloserr( int *errtype, char *losdata, int cellnum, double dlos, FILE *errfp);
#endif

