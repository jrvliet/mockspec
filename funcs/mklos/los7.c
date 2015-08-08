
// The C version of los7.f

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "files7.h"

#define NMAX 5000

int main(int argc, char *argv[]){

    // general variables
    int i, j, klos, error;
    int goodTCcells, badICcells, badDcells;
    double rdum;        
    char qsolist[80], paramlist[80], tranilist[80];
    char losdata[80], rootname[80], lostag[80], losdatafile[80];
    char new_line[200];    
    char *p;

    // los and galaxy data
    double zlos, dlos, impact;
    double Slos, vlos, vabs;
    double Rgalx, Rgaly, Rgalz, Rgal;
    double Vgalx, Vgaly, Vgalz, Vgal;
    
    // atomic data
    double mamu;
    char ionlabel[80], linesfile[80];
    
    // cell data
    int cellnum[NMAX];
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
    
    
    strcpy(qsolist, "qso.list");

    // Open the log files
    FILE *runfp = fopen("Mockspec.runlog.los7", "w");
    FILE *errfp = fopen("Mockspec.errlog.los7", "w");

    
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
    while(fgets(new_line,sizeof(new_line),listfp)){
        klos++;
    
        // Remove the return character
        p = strchr(new_line, '\n');
        *p = '\0';
        strcpy(losdata, new_line);
        
    } 


    return 0;
}


