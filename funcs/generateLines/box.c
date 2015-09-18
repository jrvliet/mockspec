
// Functions to support genLOS.c relating to the box

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "box.h"
#include "matrix.h"

int in_box(double x, double y, double z, double boxsize){
    
    int result = 0;
    double size, negsize;
    size = boxsize;
    negsize = -1.0*size;

    if(x<size && x>negsize){
        if(y<size && y>negsize){
            if(z<size && z>negsize){
                result = 1;
            }
        }
    }

    return result;
}


void find_ends(double px, double py, double pz, double dx, double dy, double dz,
    double boxsize, double *xen, double *yen, double *zen, double *xex, 
    double *yex, double *zex, double *ten, double *tex){

    int i;
    int tsize = 10000;
    int tstart = -1000;
    int tend = 1000;
    double tstep;
    double size = boxsize;
    double t[tsize];
    double x1, y1, z1, x2, y2, z2;
    
    // Intialize entry and exit points
    *xen = 0.0;
    *yen = 0.0;
    *zen = 0.0;
    *ten = 0.0;
    *xex = 0.0;
    *yex = 0.0;
    *zex = 0.0;
    *tex = 0.0;

    tstep = ((double)tend - (double)tstart)/(double)tsize;

    for(i=0; i<tsize; i++){
        t[i] = tstart + i*tstep;
    }

    for(i=0; i<tsize-1; i++){
        x1 = px + dx*t[i];
        y1 = py + dy*t[i];
        z1 = pz + dz*t[i];
        
        x2 = px + dx*t[i+1];
        y2 = py + dy*t[i+1];
        z2 = pz + dz*t[i+1];
        
        if (in_box(x1,y1,z1,size)==0 && in_box(x2,y2,z2,size)==1){
            *xen = (x1+x2)/2.0;
            *yen = (y1+y2)/2.0;
            *zen = (z1+z2)/2.0;
            *ten = (t[i]+t[i+1])/2.0;
        }
        
        if (in_box(x1,y1,z1,size)==1 && in_box(x2,y2,z2,size)==0){
            *xex = (x1+x2)/2.0;
            *yex = (y1+y2)/2.0;
            *zex = (z1+z2)/2.0;
            *tex = (t[i]+t[i+1])/2.0;
        }
    }
}


double inclination(double **a_gtb, double *db){

    // Calculates the galaxy's inclination (radians) relative to the LOS
    // Passed in: 
    //    a_gtb = rotation matrix to convert from galaxy frame to box frame    
    //    db = directional vector of LOS, defined in box frame

    // normG is a vector along the galaxy's z axis, which is along
    // the angular momentum vector. It is defined the galaxy's 
    // reference frame
    double normG[3];
    normG[0] = 0.0;
    normG[1] = 0.0;
    normG[2] = 1.0;

    // normB is normG converted to the box's frame, which the LOS is
    // defined in
    double normB[3];   
    multiply(a_gtb, normG, normB);
    
    // losLen is the length of the LOS vector 
    double losLen = mat_len(db, 3);
    
    // normLen is the legnth of normB
    double normLen = mat_len(normB, 3);

    double dotProd = dot_product(db, normB, 3);

    // Interir angle between the galaxy's angular momentum vector
    // and the LOS directional vector is given by:
    // cos(theta) = dot(a,b) / (len(a)*len(b))
    double cosTheta = dotProd / (losLen*normLen);
    double theta = acos(cosTheta);

    // Rotate so the angle is between 0 and pi/2
    while(theta>M_PI/2.0){
        theta -= M_PI/2.0;
    }

    return theta;
}    



int equals(double a, double b){

    // Determines if two floating point numbers are approximately equal
    // This is required as ANA sometimes outpus an expansion parameter
    // of 0.800 as 0.801
    
    int result = 0;
    double tol = 0.002;
    if (fabs(a-b)<=tol){
        result = 1;
    }
    
    return result;
}



void sort_LOS(double *a, int n){

    // Perfmorms a merge sort on the LOS impact parameters
    // list, a, of length n

    if (n<2) {
        return;
    }
    // Find the halfway point in the list
    int m = n/2;

    sort_LOS(a,m);
    sort_LOS(a+m, n-m);
    merge(a, n, m);
}



void merge(double *a, int n, int m){
    int i, j, k;
    double *x = calloc(n, sizeof(double));
    for (i=0, j=m, k=0; k<n; k++){
        x[k] = j==n      ? a[i++]
             : i==m      ? a[j++]
             : a[j]<a[i] ? a[j++]
             :             a[i++];
    }
    for (i=0; i<n; i++){
        a[i] = x[i];
    }
    free(x);
}




