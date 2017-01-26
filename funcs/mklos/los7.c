
// The C version of los7.f

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "files7.h"
#include "cellpathsub7.h"
#include "galframe7.h"

#define NMAX 5000

int main(int argc, char *argv[]){

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
    int losnum;
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

    // Loop over the lines of sight
    while(fgets(new_line,sizeof(new_line),listfp)){
        klos++;
        //printf("LOS Num: %d\t",klos); 
        // Remove the return character
        p = strchr(new_line, '\n');
        *p = '\0';
        strcpy(losdata, new_line);

        // Parse the file name
        cutfname(losdata, &galID, &ion, &lostag, losdatafile);
        sscanf(lostag+3, "%d", &losnum);

        mamu = getamu(tranilist, ion);

        mkfname(losdata, galID, ion, lostag, linesfile);

        // Read in the cell data for this losdata
        gethdr( klos, &a, &xg, &yg, &zg, &vxg, &vyg, &vzg, &b1, &b2, &x0, &y0, 
                &z0, &l, &m, &n, ap, &zbox, &vgal, &zgal, losdata);

        ndata = readcells( cellnum, x, y, z, vx, vy, vz, Lcell, ndencell, fion, 
                           temp, zmfrac, losdata);

        // Open the .losdata file and write the header
        losfp = fopen(losdatafile, "w");
        fprintf(losfp, "1 \t2 \t3 \t4 \t5 \t6 \t7 \t8 \t"
                       "9 \t10 \t11 \t12 \t13 \t14 \t15 \t"
                       "16 \t17 \t18 \t19 \t20 \t21 \t22 \t23 \t"
                       "24 \t25 \t26 \t27 \t28 \n"); 
        fprintf(losfp, "Slos \tRgal \tzabs \tvlos \tvabst\tdlos "
                        "\tnion \tfion \tzmfrac \tNion \tT \t"
                        "bpar \t Vgtot \tvrp \tVtheta \tVphi \tvzp "
                        "\txp \typ \tzp \trp \ttheta \tphi "
                        "\tvrp/Vgt Vth/Vgt Vph/Vgt vzp/Vg \t "
                        "cellID\n");
        fclose(losfp);
        // Inner loop: 
        // Loop over cell data for the los, compute the los quantities
        for(i=0; i<ndata; i++){

            // compute cell los quantities; Slos [kpc], Rgal [kpc], velocities
            // [km/s] columns [1/cm2], b parameters [km/s]
            Slos = l*(x[i]-xg) + m*(y[i]-yg) + n*(z[i]-zg);
            vlos = l*vx[i] + m*vy[i] + n*vz[i];
            zlos = zbox + (1.0 + zbox)*(vlos/ckms); 
            vabs = ckms * (zlos-zgal)/(1.0+zgal);

            // compute the line of sight path length through the cell (credit to
            // Bobby Edmonds' for the getD subroutine); returns the value of dlos
            // in kpc, convert to centimeters upon return
            getD( l, m, n, x0, y0, z0, x[i], y[i], z[i], Lcell[i], &dlos, 
                  &error, errtype, losnum);
            
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
            Nline[i] = ndencell[i] * dlos;
            zline[i] = zlos;
            bline[i] = 1.0e-5 * sqrt(2.0*kboltz*temp[i]/(amu*mamu));

            // Translate the coordinate system to the galaxy center
            // Translate the velocities to the galaxy velocity
            Rgalx = x[i]-xg;
            Rgaly = y[i]-yg;
            Rgalz = z[i]-zg;
            Rgal  = sqrt(Rgalx*Rgalx + Rgaly*Rgaly + Rgalz*Rgalz);
            Vgalx  = vx[i]-vxg;
            Vgaly  = vy[i]-vyg;
            Vgalz  = vz[i]-vzg;
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
            wrtlosdata( Slos, Rgal, zline[i], vlos, vabs, dlos, ndencell[i], 
                        fion[i], zmfrac[i], Nline[i], temp[i], bline[i], Vgal, 
                        vrp, V_theta, V_phi, vzp, xp, yp, zp, rp, theta, phi, 
                        cellnum[i], losdatafile);

        }  // Ends loop over cells in a single LOS

        // Write the .lines file, one for each transition
        // These files used by specsynth to generate the spectra
        wrtlines(zgal, zline, Nline, bline, cellnum, linesfile, ndata);
        //printf("Number of new cells: %d\n", ndata);

    } 
    //printf("Close runfp\n");
    fclose(runfp);
    //printf("Close errfp\n");
    fclose(errfp);
    //printf("Return\n");
    return 0;
}


