#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <math.h>

void inverse(double **start, double **end);

void multiply(double **A, double *p, double *v);

double mat_len(double *a, int n);

double dot_product(double *a, double *b, int n);

#endif
