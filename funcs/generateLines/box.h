#ifndef BOX_H
#define BOX_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int in_box(double x, double y, double z, double boxsize);

void find_ends(double px, double py, double pz, double dx, double dy, double dz, 
  double boxsize, double *xen, double *yen, double *zen, double *xex,
  double *yex, double *zex, double *ten, double *tex, int tmax);

double inclination(double **a_gtb, double *db);

int equals(double a, double b);

void sort_LOS( double *a, int n);

void merge(double *a, int n, int m);

#endif
