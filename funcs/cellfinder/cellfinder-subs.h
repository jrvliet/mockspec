#ifndef FBCELLFINDER_SUBS_H
#define FBCELLFINDER_SUBS_H

#include <stdio.h>


void read_control_file(FILE *propfp, char *gasfile, char *galID, char *rootname, double *aexpn, char *summaryLoc);

void read_summary(char *galID, double *aexpn, char *summaryLoc, double *mvir, double *rvir, double *a11, double *a12, double *a13, double *a21, double *a22, double *a23, double *a31, double *a32, double *a33, double *vpec_x, double *vpec_y, double *vpec_z);

const char * filegen(char *rootname, int losnum);
//void filegen(char *rootname, int losnum, char *outfile);

void write_OutfileHdr(FILE *outfp, double *aexpn, double *R0, double *phi, double *l, double *b, double *xen, double *yen, double *zen, double *losx, double *losy, double *losz, double *a11, double *a12, double *a13, double *a21, double *a22, double *a23, double *a31, double *a32, double *a33, double *Xcom, double *Ycom, double *Zcom, double *VXcom, double *VYcom, double *VZcom, double *x0, double *y0, double *z0, double *vx_obs, double *vy_obs,double *vz_obs);

void write_LOSprops(FILE *propsfp, int losnum, double *aexpn, double *R0, double *phi, double *l, double *b, double *xen, double *yen, double *zen, double *losx, double *losy, double *losz, double *a11, double *a12, double *a13, double *a21, double *a22, double *a23, double *a31, double *a32, double *a33, double *Xcom, double *Ycom, double *Zcom, double *VXcom, double *VYcom, double *VZcom, double *x0, double *y0, double *z0, double *vx_obs, double *vy_obs,double *vz_obs);

#endif
