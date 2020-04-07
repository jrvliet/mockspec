
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "matrix.h"
#include "box.h"

int main (int argc, char *argv[]){

    printf("Starting genLOS\n");
    int i, nLOS, ncores;
    double incDeg, maximpact, a; 
    
    // Read in command line arguments
    // 1 = galID
    // 2 = current location where code is being run
    // 3 = expansion parameter
    // 4 = inclination 
    // 5 = number of LOS
    // 6 = maximum impact paramter (in Rvir)
    // 7 = number of cores to use (not used anywhere here)
//    char *galID = argv[1];
    char *summaryLoc = argv[2];
    char *expn = argv[3];
    sscanf(expn, "%lf", &a);
    sscanf(argv[4], "%lf", &incDeg);
    sscanf(argv[5], "%d", &nLOS);
    sscanf(argv[6], "%lf", &maximpact);
    sscanf(argv[7], "%d", &ncores);

    // Convert inclination to radians
    double inc = incDeg * (M_PI/180.0);  

    int tmax = 1000;
    double b, phi, phiDeg;
    double aexpn, redshift;
    double mvir, rvir;
    double a11, a12, a13;
    double a21, a22, a23;
    double a31, a32, a33;
    double vpec_x, vpec_y, vpec_z;
    double galInc;

    double xen, yen, zen, ten;
    double xex, yex, zex, tex;
    double *xens, *yens, *zens;
    double *xexs, *yexs, *zexs;
    double db[3], pb[3], dg[3], pg[3];

    double *ds = (double *)calloc(3, sizeof(double));
    double *ps = (double *)calloc(3, sizeof(double));

    xens = (double *)calloc(nLOS, sizeof(double));
    yens = (double *)calloc(nLOS, sizeof(double));
    zens = (double *)calloc(nLOS, sizeof(double));

    xexs = (double *)calloc(nLOS, sizeof(double));
    yexs = (double *)calloc(nLOS, sizeof(double));
    zexs = (double *)calloc(nLOS, sizeof(double));

    char sumfile[300];
    strcpy(sumfile, summaryLoc);
    strcat(sumfile, "/rotmat_a");
    strcat(sumfile, expn);
    strcat(sumfile, ".txt");
    char new_line[1000];    
    FILE *fp = fopen(sumfile, "r");

    printf("summaryLoc: %s\n",summaryLoc);
    printf("expng: %s\n",expn);
    printf("sumfile: %s\n",sumfile);

    if (fp==NULL){
        printf("\nTest ERROR in genLOS\n");
        printf("Could not open %s\n", sumfile);
        printf("Exitting....\n\n");
        exit(1);
    }
    // Read in the summary file, skipping past the one header row
    fgets(new_line, sizeof(new_line), fp);
    int found = 0;
    while(fgets(new_line, sizeof(new_line), fp)){
        sscanf(new_line, "%lf %lf %lf %lf %lf "
                     "%lf %lf %lf %lf %lf "
                     "%lf %lf %lf %lf %lf "
                     "%lf", 
                      &aexpn, &redshift, &mvir, &rvir, &a11, 
                      &a12, &a13, &a21, &a22, &a23, 
                      &a31, &a32, &a33, &vpec_x, &vpec_y,
                      &vpec_z);
        if ( equals(aexpn,a)){
            found = 1;
        }
    }
    if (found==0){
        printf("\nERROR in genLOS\n");
        printf("Could not find %s in %s\n", expn, sumfile);
        printf("Exitting...\n\n");
        exit(1);
    }
    
    fclose(fp);

    

    // Seed random
    //int seed = 25525;
    //srand(seed);
    srand(time(NULL));

    double maximpact_kpc = maximpact*rvir;
    double boxsize = 4.0*rvir;

    // Open file operators
    FILE *fpdat = fopen("lines.dat", "w");
    FILE *fpinfo = fopen("lines.info", "w");
    
    // Write header
    //fprintf(fpdat, "#        Enter points               Exit Points\n");
    fprintf(fpdat, "# xen      yen         zen      xex      yex       zex\n");
    //fprintf(fpinfo, "# More details on each LOS\n");
    fprintf(fpinfo, "# LOS num        b(kpc)       phi      Inclination\n");



    // Rotation matrix to go from box to galaxy
    double **a_btg, **a_gtb, **a_stg;
    a_btg = (double **)calloc(3, sizeof(double *));
    a_gtb = (double **)calloc(3, sizeof(double *));
    a_stg = (double **)calloc(3, sizeof(double *));
    for (i=0; i<3; i++){
        a_btg[i] = (double *)calloc(3, sizeof(double));
        a_gtb[i] = (double *)calloc(3, sizeof(double));
        a_stg[i] = (double *)calloc(3, sizeof(double));
    }
        
    // Fill a_btg, the rotation matrix for translating from
    // the box frame to the galaxy frome. This is the rotation
    // matrix output by ANA and stored in the summary file
    a_btg[0][0] = a11;
    a_btg[0][1] = a12;
    a_btg[0][2] = a13;
    
    a_btg[1][0] = a21;
    a_btg[1][1] = a22;
    a_btg[1][2] = a23;

    a_btg[2][0] = a31;
    a_btg[2][1] = a32;
    a_btg[2][2] = a33;
    
    // Get the galaxy to box rotation matrix, which is the
    // inversion of the box to galaxy rotation matrix
    inverse(a_btg, a_gtb);

    // Fill ds, the directional vector for the LOS in the sky fram
    ds[0] = 0.0;
    ds[1] = 0.0;
    ds[2] = -1.0;

    // Fill a_stg, the rotation matrix to go from sky to galaxy
    // Depends only on inclination
    a_stg[0][0] = 1.0;    
    a_stg[0][1] = 0.0;    
    a_stg[0][2] = 0.0;    

    a_stg[1][0] = 0.0;    
    a_stg[1][1] = cos(inc);    
    a_stg[1][2] = -1.0*sin(inc);    

    a_stg[2][0] = 0.0;    
    a_stg[2][1] = sin(inc);    
    a_stg[2][2] = cos(inc);    

    // Find dg, the direction vecotr for the LOS in the galaxy frame
    // This is a_stg*ds
    multiply(a_stg, ds, dg);

    // Generate impact parameters
    double impacts[nLOS];
    for(i=0; i<nLOS; i++){
        b = ((double)rand()/(double)RAND_MAX)*maximpact_kpc;
        impacts[i] = b;
    }
    // Sort the impact paramters
    sort_LOS(impacts, nLOS);

    // Loop over LOS
    for(i=0; i<nLOS; i++){

        b = impacts[i];
        phi = ((double)rand()/(double)RAND_MAX)*(2.0*M_PI);
        phiDeg = phi * (180.0/M_PI);
        fprintf(fpinfo, "%d \t %lf \t %lf \t %lf\n",i+1,b,phiDeg,incDeg);

        // Fill ps, the position vector in the sky frame
        ps[0] = b*cos(phi);
        ps[1] = b*sin(phi);
        ps[2] = 0.0;

        // Translate the position vector to the galaxy frame
        // pg = a_stg*ps
        multiply(a_stg, ps, pg);
        
        // Rotate the galaxy's directional vector and impact point 
        // into the box's frame
        // db = a_gtb*dg
        // pb = a_gtb*pg
        multiply(a_gtb, dg, db);
        multiply(a_gtb, pg, pb);
        
        // Find the enter and exit point of the LOS 
        find_ends(pb[0], pb[1], pb[2], db[0], db[1], db[2], boxsize,
                  &xen, &yen, &zen, &xex, &yex, &zex, &ten, &tex, tmax);

        
        if ((xex==0.) || (yex==0.) || (zex==0.) || 
           (xen==0.) || (yen==0.) || (zen==0.)){
            tmax *= 20;
            find_ends(pb[0], pb[1], pb[2], db[0], db[1], db[2], boxsize,
                      &xen, &yen, &zen, &xex, &yex, &zex, &ten, &tex, tmax);
        }


        xens[i] = xen;
        yens[i] = yen;
        zens[i] = zen;
        
        xexs[i] = xex;
        yexs[i] = yex;
        zexs[i] = zex;

        fprintf(fpdat, "%lf \t%lf \t%lf \t%lf \t%lf \t%lf\n",xen,yen,zen,xex,yex,zex);
//        printf("%lf \t%lf \t%lf \t%lf \t%lf \t%lf\n",xen,yen,zen,xex,yex,zex);
//        printf("\n");
    }

    // Test to ensure the inclination of the galaxy is actually what 
    // inc is set to
    galInc = inclination(a_gtb, db);     //radians
    // Report the inclination to screen
//    printf("Calculated inclination angle: %lf degrees\n", galInc*(180.0/M_PI));


    FILE *fptest = fopen("testing.dat", "a");
    fprintf(fptest, "%lf \t %lf\n", inc, galInc);
    fclose(fptest);

    // Close all files
    fclose(fpdat);
    fclose(fpinfo);

    return 0;
}







