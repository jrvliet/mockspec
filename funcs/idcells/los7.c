
// The C version of los7.f
// This version is a functional form that is called by idcells.c
// It operates one cell at a time

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "files7.h"
#include "cellpathsub7.h"
#include "galframe7.h"

#define NMAX 5000



int los7(char *galID, char *ion, char *losname, double mamu, struct los losprops, 
            struct cell cprops, struct orient oprops, stuct gal gprops){

    // general variables
    int i, klos, error, errtype[3];
    int badDcells;

    char *qsolist = calloc(80, sizeof(char));
    char *paramlist = calloc(80, sizeof(char));
    char *tranilist = calloc(80, sizeof(char));
    char *ion = calloc(80, sizeof(char));
    char *losdata = calloc(80, sizeof(char));
    char *galID = calloc(80, sizeof(char));
    char *lostag = calloc(80, sizeof(char));
    char *losdatafile = calloc(80, sizeof(char));

    char new_line[200];    
    char *p;
    FILE *losfp;

    // Constants
    double ckms = 2.99792458e5;
    double pc2cm = 3.261633*9.460528e17;
    double kpc2cm = 1000.0 * pc2cm;
    double amu     = 1.66053878e-24;
    double kboltz  = 1.380658e-16;

    // los and galaxy data
    double zlos, dlos; //, impact;
    double Slos, vlos, vabs;
    double Rgalx, Rgaly, Rgalz, Rgal;
    double Vgalx, Vgaly, Vgalz, Vgal;
    
    // atomic data
    double mamu;
    char linesfile[80];
    
    // cell data
    int cellnum[NMAX], ndata;
    double x[NMAX], y[NMAX], z[NMAX];
    double vx[NMAX], vy[NMAX], vz[NMAX];
    double ndencell[NMAX], fion[NMAX], zmfrac[NMAX];
    double temp[NMAX], Lcell[NMAX];
    double Nline[NMAX], zline[NMAX], bline[NMAX];
    
    // box reference frame quantities
    double a, zbox, x0, y0, z0;
    double l, m, n, b1, b2;
    double xg, yg, zg, vxg, vyg, vzg;
    double vgal, zgal;
    
    // galaxy reference frame quantities
    double ap[3][3];
    double xp, yp, zp, rp, theta, phi;
    double vxp, vyp, vzp, vrp, V_theta, V_phi;
    
    
    // Set up runtime files
    strcpy(tranilist, "Mockspec.transitions");
    strcpy(paramlist, "Mockspec.runpars");
    strcpy(qsolist, "qso.list");

    // Open the log files
    // FILE *runfp = fopen("Mockspec.runlog.los7", "w");
    FILE *errfp = fopen("Mockspec.errlog.los7", "a");

    
    // Open the list file of LOS files output by ana.2 (Daniel's code)
    // loop over this list file; for each LOS, raw cell quantities are
    // stored in arrays, whereas we will compute LOS data for each cell
    // on the fly
    FILE *listfp = fopen(qsolist, "r");
    if (listfp==NULL){
        printf("\n\nERROR in los7.c\n");
        printf("Cannot open file %s\n", qsolist);
        printf("Exiting...\n");
        exit(1);
    }

    klos = 0; 

    x = cprops.x;
    y = cprops.y;
    z = cprops.z;
    vx = cprops.vx;
    vy = cprops.vy;
    vz = cprops.vz;
    Lcell = cprops.size;
    cellnum = cprops.num;

    x0 = lprops.x0;
    y0 = lprops.y0;
    z0 = lprops.z0;
    
    vxg = gprops.vpec_x;
    vyg = gprops.vpec_y;
    vzg = gprops.vpec_z;

    l = oprops.losx;
    m = oprops.loxy;
    n = oprops.losz;
    xg = 0.0;
    yg = 0.0;
    zg = 0.0;
   
    ndencell = cprops.nden;
    temp = cprops.temp; 
    
    // compute cell los quantities; Slos [kpc], Rgal [kpc], velocities
    // [km/s] columns [1/cm2], b parameters [km/s]
    Slos = l*(x-xg) + m*(y-yg) + n*(z-zg);
    vlos = l*vx + m*vy + n*vz;
    zlos = zbox + (1.0 + zbox)*(vlos/ckms); 
    vabs = ckms * (zlos-zgal)/(1.0+zgal);

    // compute the line of sight path length through the cell (credit to
    // Bobby Edmonds' for the getD subroutine); returns the value of dlos
    // in kpc, convert to centimeters upon return
    getD( l, m, n, x0, y0, z0, x[i], y[i], z[i], Lcell[i], &dlos, 
          &error, errtype);
    
    // If getD returns clean, convert from Mpc to kpc
    // If getD returns an error, use the cube root of the
    // cell volume; communicate to the error log file
    if (error==1){
        dlos = kpc2cm*Lcell[i];
        badDcells = badDcells+1;
        dloserr(errtype, losdata, cellnum[i], dlos, errfp);    
    }
    else{
        dlos = kpc2cm*dlos;
    }

    // Compute the cell column density [cm^-2]
    // observed redshift
    // Doppler param [km/s]
    Nline = ndencell * dlos;
    zline = zlos;
    bline = 1.0e-5 * sqrt(2.0*kboltz*temp/(amu*mamu));

    // Translate the coordinate system to the galaxy center
    // Translate the velocities to the galaxy velocity
    Rgalx = x-xg;
    Rgaly = y-yg;
    Rgalz = z-zg;
    Rgal  = sqrt(Rgalx*Rgalx + Rgaly*Rgaly + Rgalz*Rgalz);
    Vgalx  = vx-vxg;
    Vgaly  = vy-vyg;
    Vgalz  = vz-vzg;
    Vgal   = sqrt(Vgalx*Vgalx + Vgaly*Vgaly + Vgalz*Vgalz);

    // Rotate physical coordinates (spatial will be kpc)
    // Rotate velocity vectors
    // Compute spherical coordinate velocities
    rp = 0.0;
    rotate(ap, Rgalx, Rgaly, Rgalz, &xp, &yp, &zp);
    rotate(ap, Vgalx, Vgaly, Vgalz, &vxp, &vyp, &vzp);
    sphvels(xp, yp, zp, &rp, &theta, &phi, 
            vxp, vyp, vzp, &vrp, &V_theta, &V_phi);

    // Write processed data to the .losdata file
    wrtlosdata( Slos, Rgal, zline, vlos, vabs, dlos, ndencell, 
                fion, zmfrac, Nline, temp, bline, Vgal, 
                vrp, V_theta, V_phi, vzp, xp, yp, zp, rp, theta, phi, 
                cellnum, losdatafile);

    // Write the .lines file, one for each transition
    // These files used by specsynth to generate the spectra
    wrtlines(zgal, zline, Nline, bline, cellnum, linesfile, ndata);


    return 0;
}


