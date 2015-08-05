#include "lmn_cellfinder.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void lmn(){
  printf("\n");
}

void POBox(){
  printf("\n");
}

void paradox_chk(){
  printf("\n");
}

void get_los_endpt(){
  printf("\n");
}

void lmn4(double xen, double yen, double zen, double xex, double yex, double zex, double Xcom, double Ycom, double Zcom, double *losx, double *losy, double *losz, double *x0, double *y0, double *z0){

  double dx, dy, dz, denom;
  double lb, mb, nb, dl;
  
  // Define the directional vector
  dx = xex-xen;
  dy = yex-yen;
  dz = zex-zen;

  // Determine the directional cosines
  denom = sqrt(dx*dx + dy*dy + dz*dz);
  lb = dx / denom;
  mb = dy / denom;
  nb = dz / denom;

  // Determine a unit vector describing the LOS
  //# Normalize the directional vector
  dl = sqrt(dx*dx + dy*dy + dz*dz);
  dx = dx/dl;
  dy = dy/dl;
  dz = dz/dl;


  // Parametric version with t=1
  *x0 = xen + dx;
  *y0 = yen + dy;
  *z0 = zen + dz;
  
  // Return the directional cosines
  *losx = lb;
  *losy = mb;  
  *losz = nb;
  
  

}
