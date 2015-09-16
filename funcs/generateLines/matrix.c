
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"



void inverse(double **start, double **end){
    // Finds the inverse of a 3x3 matrix

    double a, b, c, d, e, f, g, h, i;
    double A, B, C, D, E, F, G, H, I;
    double det;

    a = start[0][0];
    b = start[0][1];
    c = start[0][2];
    
    d = start[1][0];
    e = start[1][1];
    f = start[1][2];
    
    g = start[2][0];
    h = start[2][1];
    i = start[2][2];

    // Find the Cofactors
    A = e*i - f*h;
    B = -1.0 * (d*i - f*g);
    C = d*h - e*g;
    D = -1.0*(b*i - b*h);
    E = a*i - c*g;
    F = -1.0*(a*h - b*g);
    G = b*f - c*e;
    H = -1.0*(a*f - c*d);
    I = a*e - b*d;
    
    // Get the determinate from the rule of Sarrus
    det = a*A + b*B + c*C;
    
    // Set the ending matrix
    end[0][0] = A/det;
    end[0][1] = D/det;
    end[0][2] = G/det;
    
    end[1][0] = B/det;
    end[1][1] = E/det;
    end[1][2] = H/det;
    
    end[2][0] = C/det;
    end[2][1] = F/det;
    end[2][2] = I/det;
    
}



void multiply(double **A, double *p, double *v){

    // Mulitplies A*p, outputs result to v
    // A = 3x3
    // p = 3x1
    // v = 3x1

    // A = [ a b c ]     p = [ x ]     
    //     [ d e f ]         [ y ]
    //     [ g h i ]         [ z ]

    double a, b, c, d, e, f, g, h, i, x, y, z;

    a = A[0][0];
    b = A[0][1];
    c = A[0][2];
    
    d = A[1][0];
    e = A[1][1];
    f = A[1][2];
    
    g = A[2][0];
    h = A[2][1];
    i = A[2][2];

    x = p[0];
    y = p[1];
    z = p[2];

    v[0] = a*x + b*y + c*z;
    v[1] = d*x + e*y + f*z;
    v[2] = g*x + h*y + i*z;

}


double mat_len(double *a, int n){

    // Computes the length of the matrix a with n entris

    double len = 0.0;
    int i;
    
    for (i=0; i<n; i++){
        len += a[i]*a[i];
    }
    len = sqrt(len);

    return len;
}



double dot_product(double *a, double *b, int n){

    // Computes the dot product of a and b, which are 
    // matrixes of length n

    double result = 0.0;
    int i;

    for (i=0; i<n; i++){
        result += a[i]*b[i];
    }
    
    return result;
}


        





