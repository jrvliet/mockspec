#ifndef GALFRAME7_H
#define GALFRAME7_H

#include <stdio.h>

void rotate(double a[][3], double x, double y, double z, 
            double *xp, double *yp, double *zp);

void sphvels( double xp, double yp, double zp, double *rp, double *theta, 
              double *phi, double vxp, double vyp, double vzp, double *vrp, 
              double *vtheta, double *vphi);

#endif
