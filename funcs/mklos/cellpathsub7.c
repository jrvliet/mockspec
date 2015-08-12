

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//     DESCRIPTION:
//     this subroutine determines the actual pathlength that the los
//     pierces through the cell; FYI it is not an efficient algorithm,
//     but all tests indicate that it works; there is a condition for
//     failure and errors are returned in such cases

void printArr( double *x, int size){
    int i;
    for (i=0; i<size; i++){
        printf("\t%lf", x[i]);
    }
    printf("\n");
}


void getD(double l, double m, double n, double x0, double y0, double z0, 
          double X, double Y, double Z, double lengthcell, double *cellpath, 
          int *error, int *errtype){

    double lq, mq, nq;
    double xu[4], xd[4], yu[4], yd[4], zu[4], zd[4];
    double pnt1[4], pnt2[4];//, pnt3[4]
    double diffpnt[4];
    double a, k;

    // Initialize 
    *error = 0;
    errtype[0] = 0;
    errtype[1] = 0;
    errtype[0] = 0;


    // The parametric distance along the LOS 
    // (What everyone else in world calls t)
    a = 0.0;

    // Various flags for the directoin cosines, to check if the
    // LOS is perpendicular to any of the cell walls
    lq = 0.0;
    mq = 0.0;
    nq = 0.0;
    k = 0.0;    // Not used in practice

    // These will be used to compute cellpath quantity 
    // LOS entry points on cell wall
    pnt1[0] = 0.0;
    pnt1[1] = 0.0;
    pnt1[2] = 0.0;
    pnt1[3] = 0.0;

    // LOS exit points on cell wall
    pnt2[0] = 0.0;
    pnt2[1] = 0.0;
    pnt2[2] = 0.0;
    pnt2[3] = 0.0;

    // Difference between entry and exit points
    diffpnt[0] = 0.0;
    diffpnt[1] = 0.0;
    diffpnt[2] = 0.0;
    diffpnt[3] = 0.0;   

    *cellpath = 0.0;
    
//    printf("\nIn getD:\n");
    
//    printf("\tl: %lf \n\tm: %lf \n\t n: %lf \n", l, m, n);
    
    // Sets flags to check if LOS is perpidicular to a cell wall
    // if not perpidicular then set flag high
    if (l>0.0){
        lq = 1.0;
    }
    if (l<0.0){
        lq = 1.0;
    }
    
    if (m>0.0){
        mq = 1.0;
    }
    if (m<0.0){
        mq = 1.0;
    }

    if (n>0.0){
        nq = 1.0;
    }
    if (n<0.0){
        nq = 1.0;
    }

//    printf("\tlq: %lf \n\tmq: %lf \n\t nq: %lf \n", lq, mq, nq);

    // Following code initializes points of intersection on cell
    // Initializing the matrix of points on each face of the cell

    xu[0] = 0.0;    // x center x pos in front of cell center
    xu[1] = 0.0;    // x center y pos in front of cell center
    xu[2] = 0.0;    // x center z pos in front of cell center
    xu[3] = 1.0;    // flag x center in front of cell center

    xd[0] = 0.0;    // x center x pos behind cell center
    xd[1] = 0.0;    // x center y pos behind cell center
    xd[2] = 0.0;    // x center z pos behind cell center
    xd[3] = 1.0;    // flag x center behind cell center

    yu[0] = 0.0;    // y center x pos in front of cell center
    yu[1] = 0.0;    // y center y pos in front of cell center
    yu[2] = 0.0;    // y center z pos in front of cell center
    yu[3] = 1.0;    // flag y center in front of cell center

    yd[0] = 0.0;    // y center x pos behind cell center
    yd[1] = 0.0;    // y center y pos behind cell center
    yd[2] = 0.0;    // y center z pos behind cell center
    yd[3] = 1.0;    // flag y center behind cell center

    zu[0] = 0.0;    // z center x pos in front of cell center
    zu[1] = 0.0;    // z center y pos in front of cell center
    zu[2] = 0.0;    // z center z pos in front of cell center
    zu[3] = 1.0;    // flag z center in front of cell center

    zd[0] = 0.0;    // z center x pos behind cell center
    zd[1] = 0.0;    // z center y pos behind cell center
    zd[2] = 0.0;    // z center z pos behind cell center
    zd[3] = 1.0;    // flag z center behind cell center

    
    // following code finds possible points of intersection
    // draw line from x0,y0,z0 position of line of sight (on BOX)
    // to find the xu(1,2,3), xd(1,2,3) points on the cell where
    // it is intersected by the line of sight

    // if not perpendicular then
    // xu = the probable x center x,y,z pos in front of the cell
    // xd = the probable x center x,y,z pos behind the cell
    if (lq>0.0){
        xu[0] = X + (lengthcell/2.0);
        a = (xu[0]-x0)/l;        // distance along LOS from X0, Y0, Z0
        xu[1] = y0 + m*a;
        xu[2] = z0 + n*a;

        xd[0] = X - (lengthcell/2.0);
        a = (xd[0]-x0) / l;
        xd[1] = y0 + m*a;
        xd[2] = z0 + n*a;
    }

    // if not perpendicular then
    // yu = the probable y center x,y,z pos in front of the cell
    // yd = the probable y center x,y,z pos behind the cell
    if (mq > 0.0) {
        yu[1] = Y + (lengthcell/2.0);
        a = (yu[1]-y0)/m;
        yu[0] = x0 + l*a;
        yu[2] = z0 + n*a;
        yd[1] = Y - (lengthcell/2.0);
        a = (yd[1]-y0)/m;
        yd[0] = x0 + l*a;
        yd[2] = z0 + n*a;
    }
    
    // if not perpendicular then
    // yu = the probable y center x,y,z pos in front of the cell
    // yd = the probable y center x,y,z pos behind the cell
    if (nq>0.0) {
        zu[2] = Z + (lengthcell/2.0);
        a = (zu[2]-z0)/n;
        zu[0] = x0 + l*a;
        zu[1] = y0 + m*a;
        zd[2] = Z - (lengthcell/2.0);
        a = (zd[2]-z0)/n;
        zd[0] = x0 + l*a;
        zd[1] = y0 + m*a;
    }
   
    
    // Following code test whether points are on cell surface
    // If outside the cell, then set flag low
    // Check both fornt and back cell walls
    
//    printf("\nX: %lf\tY: %lf\tZ: %lf\n\n", X, Y, Z);
    // X posiitions
    if( (xu[1] > (Y+(lengthcell/2.0))) || (xu[1] < (Y-(lengthcell/2.0))) ) {
        xu[3] = 0.0;
    }
    if( (xu[2] > (Z+(lengthcell/2.0))) || (xu[2] < (Z-(lengthcell/2.0))) ) {
        xu[3] = 0.0;
    }
    if( (xd[1] > (Y+(lengthcell/2.0))) || (xd[1] < (Y-(lengthcell/2.0))) ) {
        xd[3] = 0.0;
    }
    if( (xd[2] > (Z+(lengthcell/2.0))) || (xd[2] < (Z-(lengthcell/2.0))) ) {
        xd[3] = 0.0;
    }
      

    // Y posiitions
    if( (xu[0] > (X+(lengthcell/2.0))) || (yu[0] < (X-(lengthcell/2.0))) ) {
        yu[3] = 0.0;
    }
    if( (yu[2] > (Z+(lengthcell/2.0))) || (yu[2] < (Z-(lengthcell/2.0))) ) {
        yu[3] = 0.0;
    }
    if( (yd[0] > (X+(lengthcell/2.0))) || (yd[0] < (X-(lengthcell/2.0))) ) {
        yd[3] = 0.0;
    }
    if( (yd[2] > (Z+(lengthcell/2.0))) || (yd[2] < (Z-(lengthcell/2.0))) ) {
        yd[3] = 0.0;
    }
      

    // Z posiitions
    if( (zu[1] > (Y+(lengthcell/2.0))) || (zu[1] < (Y-(lengthcell/2.0))) ) {
        zu[3] = 0.0;
    }
    if( (zu[0] > (X+(lengthcell/2.0))) || (zu[0] < (X-(lengthcell/2.0))) ) {
        zu[3] = 0.0;
    }
    if( (zd[1] > (Y+(lengthcell/2.0))) || (zd[1] < (Y-(lengthcell/2.0))) ) {
        zd[3] = 0.0;
    }
    if( (zd[0] > (X+(lengthcell/2.0))) || (zd[0] < (X-(lengthcell/2.0))) ) {
        zd[3] = 0.0;
    }
      

    // Following code extracts the entry points that are on cell surface
 /*   printf("\nxu\n");
    printArr(xu, 4);
    printf("yu\n");
    printArr(yu, 4);
    printf("zu\n");
    printArr(zu, 4);

    printf("xd\n");
    printArr(xu, 4);
    printf("yd\n");
    printArr(yd, 4);
    printf("zd\n");
    printArr(zd, 4);
    printf("\n");   */
    if (xu[3]>0.5){             // If on cell wall
        if (k<0.5){             // Appears to have no use
            pnt1[0] = xu[0];    // Store the x pos
            pnt1[1] = xu[1];    // Store the y pos
            pnt1[2] = xu[2];    // Store the z pos
            pnt1[3] = xu[3];    // Store the current value of the flag
            k = 1.0;
            xu[3] = 0.0;        // Reset flag low for additional checks below
        }
    }
    if (xd[3]>0.5){             // If on cell wall
        if (k<0.5){             // Appears to have no use
            pnt1[0] = xd[0];    // Store the x pos
            pnt1[1] = xd[1];    // Store the y pos
            pnt1[2] = xd[2];    // Store the z pos
            pnt1[3] = xd[3];    // Store the current value of the flag
            k = 1.0;
            xd[3] = 0.0;        // Reset flag low for additional checks below
        }
    }
    if (yu[3]>0.5){             // If on cell wall
        if (k<0.5){             // Appears to have no use
            pnt1[0] = yu[0];    // Store the x pos
            pnt1[1] = yu[1];    // Store the y pos
            pnt1[2] = yu[2];    // Store the z pos
            pnt1[3] = yu[3];    // Store the current value of the flag
            k = 1.0;
            yu[3] = 0.0;        // Reset flag low for additional checks below
        }
    }
    if (yd[3]>0.5){             // If on cell wall
        if (k<0.5){             // Appears to have no use
            pnt1[0] = yd[0];    // Store the x pos
            pnt1[1] = yd[1];    // Store the y pos
            pnt1[2] = yd[2];    // Store the z pos
            pnt1[3] = yd[3];    // Store the current value of the flag
            k = 1.0;
            yd[3] = 0.0;        // Reset flag low for additional checks below
        }
    }
    if (zu[3]>0.5){             // If on cell wall
        if (k<0.5){             // Appears to have no use
            pnt1[0] = zu[0];    // Store the x pos
            pnt1[1] = zu[1];    // Store the y pos
            pnt1[2] = zu[2];    // Store the z pos
            pnt1[3] = zu[3];    // Store the current value of the flag
            k = 1.0;
            zu[3] = 0.0;        // Reset flag low for additional checks below
        }
    }
    if (zd[3]>0.5){             // If on cell wall
        if (k<0.5){             // Appears to have no use
            pnt1[0] = zd[0];    // Store the x pos
            pnt1[1] = zd[1];    // Store the y pos
            pnt1[2] = zd[2];    // Store the z pos
            pnt1[3] = zd[3];    // Store the current value of the flag
            k = 1.0;
            zd[3] = 0.0;        // Reset flag low for additional checks below
        }
    }

    // Following code extracts the exit points that are on cell surface
    if (xu[3]>0.5){             // If on cell wall
        pnt2[0] = xu[0];    // Store the x pos
        pnt2[1] = xu[1];    // Store the y pos
        pnt2[2] = xu[2];    // Store the z pos
        pnt2[3] = xu[3];    // Store the current value of the flag
    }
    if (xd[3]>0.5){             // If on cell wall
        pnt2[0] = xd[0];    // Store the x pos
        pnt2[1] = xd[1];    // Store the y pos
        pnt2[2] = xd[2];    // Store the z pos
        pnt2[3] = xd[3];    // Store the current value of the flag
    }
    if (yu[3]>0.5){             // If on cell wall
        pnt2[0] = yu[0];    // Store the x pos
        pnt2[1] = yu[1];    // Store the y pos
        pnt2[2] = yu[2];    // Store the z pos
        pnt2[3] = yu[3];    // Store the current value of the flag
    }
    if (yd[3]>0.5){             // If on cell wall
        pnt2[0] = yd[0];    // Store the x pos
        pnt2[1] = yd[1];    // Store the y pos
        pnt2[2] = yd[2];    // Store the z pos
        pnt2[3] = yd[3];    // Store the current value of the flag
    }
    if (zu[3]>0.5){             // If on cell wall
        pnt2[0] = zu[0];    // Store the x pos
        pnt2[1] = zu[1];    // Store the y pos
        pnt2[2] = zu[2];    // Store the z pos
        pnt2[3] = zu[3];    // Store the current value of the flag
    }
    if (zd[3]>0.5){             // If on cell wall
        pnt2[0] = zd[0];    // Store the x pos
        pnt2[1] = zd[1];    // Store the y pos
        pnt2[2] = zd[2];    // Store the z pos
        pnt2[3] = zd[3];    // Store the current value of the flag
    }

    // Following code calculates distance between points
    // Returns 0 if basic error analysis determines if error occurred
    diffpnt[0] = pnt1[0] - pnt2[0];
    diffpnt[1] = pnt1[1] - pnt2[1];
    diffpnt[2] = pnt1[2] - pnt2[2];
    diffpnt[3] = pnt1[3] - pnt2[3];     // Why? Unknown
/*        
    printf("Point 1: \n");
    printArr(pnt1, 3);
    printf("Point 2: \n");
    printArr(pnt2, 3);
    printf("Diff Point:\n");
    printArr(diffpnt, 3); 
*/
    *cellpath = pow(diffpnt[0], 2.0) + pow(diffpnt[1], 2.0) + pow(diffpnt[2], 2.0) ;
    *cellpath = sqrt(*cellpath);

    // Final Checks
    // Check 1
    if (pnt1[0]==0.0 || pnt1[1]==0.0 || pnt1[3]==0.0){
        *cellpath = 0.0;
        errtype[0] = 1;
    }
    
    // Check 2
    if (pnt2[0]==0.0 || pnt2[1]==0.0 || pnt2[3]==0.0){
        *cellpath = 0.0;
        errtype[1] = 1;
    }

    // Check 3
    if (xu[3] + xd[3] + yu[3] + yd[3] + zu[3] + zd[3] > 1.5){
        *cellpath = 0.0;
        errtype[2] = 1;
    }

//    printf("Cellpath: %lf\n", *cellpath);
//    printf("Error: %d\n", *error);

    if (*cellpath==0.0){
        *error = 1;
    }


}
