
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Cuts up the string infile to get the galID, ion, los number, and
// create the output file name
void cutfname(char *infile, char *galID, char *ion, char *lostag, char *outfile){

    const char delim[2] = ".";
    char *input;
    // Format of infile:
    // galID.ion.los####.dat

    strcpy(input, infile);
    galID = strtok(input, delim);
    ion = strtok(input, delim);
    lostag = strtok(input, delim);      
        
    // Format of outfile
    // galID.ion.los####.losdata
    strcpy(outfile, galID);
    strcat(outfile, ".");
    strcat(outfile, ion);
    strcat(outfile, ".");
    strcat(outfile, lostag);
    strcat(outfile, ".losdata");

}



//    read in the header for each *.dat file using the parzeline
//    subroutine, which is in the module "parze.f"
//    each line of the header is read as a single string of length 200
//    characters; the string is then parzed into 'fields" separated by
//    one or more spaces.  the fields are stored as double precision
//    reals (if a field is a string, a value of 0.0d0 is returned in
//    that field).
//    if Nmaxfields needs increasing, it can be done so in the 'los7.h"
//    file
//

void gethdr( int klos, double *a, double *xg, double *yg, double *zg, double *vxg, double *vyg, double *vzg, double *b1, double *b2, double *x0, double *y0, double *z0, double *l, double *m, double *n, double ap[][3], double *zbox, double *vgal, double *zgal, char *infile){

    int i;
    double ckms = 3.0e5;
    char new_line[500];
    char *token;
    FILE *fp = fopen(infile, "r");
    double a0, xg0, yg0, zg0, vxg0, vyg0, vzg0;
    double b10, b20, Robs, phiobs;
    double x00, y00, z00, l0, m0, n0;
    double zbox0, vgal0, zgal0;     


    // First line
    fgets(new_line, sizeof(new_line), fp);
    token = strtok(new_line, " ");
    a0 = atof( strtok(new_line, " ") );     // Expansion parameter
    for (i=0;i<4;i++){
        token = strtok(new_line, " ");
    }
    xg0 = atof( strtok(new_line, " ") );    // Galaxy x center in Mpc
    yg0 = atof( strtok(new_line, " ") );    // Galaxy y center in Mpc
    zg0 = atof( strtok(new_line, " ") );    // Galaxy z center in Mpc
    for (i=0; i<5; i++){
        token = strtok(new_line, " ");
    }
    vxg0 = atof( strtok(new_line, " ") );   // Galaxy vx km/s
    vyg0 = atof( strtok(new_line, " ") );   // Galaxy vy km/s
    vzg0 = atof( strtok(new_line, " ") );   // Galaxy vz km/s

    // Second line
    fgets(new_line, sizeof(new_line), fp);
    for (i=0; i<5; i++){
        token = strtok(new_line, " ");
    }
    b10 = atoi( strtok(new_line, " ") );    // Galactic longitude
    b20 = atoi( strtok(new_line, " ") );    // Galactic latitude
    for (i=0; i<3; i++){
        token = strtok(new_line, " "); 
    }
    Robs = atof( strtok(new_line, " ") );     // R observer on Gal disk (should be 0)
    phiobs = atof( strtok(new_line, " ") );   // Azimuthal pos of obs on disk (should be 90)

    // Third line 
    // Nothing on this line
    fgets(new_line, sizeof(new_line), fp);
    
    // Forth line
    // x0, y0, z0, l, m, n
    fgets(new_line, sizeof(new_line), fp);
    for (i=0; i<6; i++){
        token = strtok(new_line, " ");
    }
    x00 = atof( strtok(new_line, " ") );    // Box entry point
    y00 = atof( strtok(new_line, " ") );
    z00 = atof( strtok(new_line, " ") );
    for (i=0; i<4; i++){
        token = strtok(new_line, " ");
    }
    l0 = atof( strtok(new_line, " ") );     // Directional Cosines
    m0 = atof( strtok(new_line, " ") );
    n0 = atof( strtok(new_line, " ") );
    
    // Fifth, sixth, seventh lines
    // Rotation matrix
    fgets(new_line, sizeof(new_line), fp);
    for (i=0; i<8; i++){
        token = strtok(new_line, " ");
    }
    ap[0][0] = atof( strtok(new_line, " ") );
    ap[0][1] = atof( strtok(new_line, " ") );
    ap[0][2] = atof( strtok(new_line, " ") );
    
    fgets(new_line, sizeof(new_line), fp);
    for (i=0; i<8; i++){
        token = strtok(new_line, " ");
    }
    ap[1][0] = atof( strtok(new_line, " ") );
    ap[1][1] = atof( strtok(new_line, " ") );
    ap[1][2] = atof( strtok(new_line, " ") );

    fgets(new_line, sizeof(new_line), fp);
    for (i=0; i<8; i++){
        token = strtok(new_line, " ");
    }
    ap[2][0] = atof( strtok(new_line, " ") );
    ap[2][1] = atof( strtok(new_line, " ") );
    ap[2][2] = atof( strtok(new_line, " ") );


    // Done reading the file header

    // Compute the redshift and the galaxy los velocity
    zbox0 = 1.0/a0 - 1.0;
    vgal0 = l0*vxg0 + m0*vyg0 + n0*vzg0;
    zgal0 = zbox0 + (1.0+zbox0)*vgal0/ckms;

    // Assign the values to the passed in pointers
    *a = a0;
    *xg = xg0;
    *yg = yg0;
    *zg = zg0;
    *vxg = vxg0;
    *vyg = vyg0;
    *vzg = vzg0;
    *b1 = b10;
    *b2 = b20;
//    *Robs = Robs0;
//    *phiobs = phiobs0;
    *x0 = x00;
    *y0 = y00;
    *z0 = z00;
    *l = l0;
    *m = m0;
    *n = n0;
    *zbox = zbox0;
    *vgal = vgal0;
    *zgal = zgal0;

    fclose(fp);
}





void readcells( int *cellnum, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *Lcell, double *ndencell, double *fion, double *temp, double *zmfrac, char *infile){

    
//     now read in the cell data and store in arrays all cell quantities
//     in proper units

    int i, cellID;
    char new_line[500];
    int hdrlen = 10;    
    double size, xLoc, yLoc, zLoc, velx, vely, velz;
    double n, t, snIa, snII, natom, fion0, nion0;
    
//     now read in the cell data and store in arrays
//     all cell quantities in proper units
//     Lcell(i)    = cell length [kpc]
//     x(i)        = cell x box position [kpc]
//     y(i)        = cell y box position [kpc]
//     z(i)        = cell z box position [kpc]
//     vx(i)       = cell x velocity component [km/s]
//     vy(i)       = cell y velocity component [km/s]
//     vy(i)       = cell z velocity component [km/s]
//     ndencell(i) = ion density [atoms/ cm^-3]
//     temp(i)     = temperature [Kelvin]
//     ZSNII       = Type II SN metal mass fraction
//     ZSNIa   ,   = Type Ia SN metal mass fraction
//     cellnum(i)  = cell number ID

    FILE *fp = fopen(infile, "r");
    
    // Loop over the header
    for (i=0; i<hdrlen; i++){
        fgets(new_line, sizeof(new_line), fp);
    }

    // Loop over all cells in the file
    i=0;
    while(fgets(new_line, sizeof(new_line), fp)){        
        
        sscanf(new_line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", &size, &xLoc, &yLoc, &zLoc, &velx, &vely, &velz, &n, &t, &snIa, &snII, &natom, &fion0, &nion0, &cellID);
        /*
        size = atof(strtok(new_line, " "));
        xLoc = atof(strtok(new_line, " "));
        yLoc = atof(strtok(new_line, " "));
        zLoc = atof(strtok(new_line, " "));
        velx = atof(strtok(new_line, " "));
        vely = atof(strtok(new_line, " "));
        velz = atof(strtok(new_line, " "));
        n = atof(strtok(new_line, " "));
        t = atof(strtok(new_line, " "));
        snIa = atof(strtok(new_line, " "));
        snII = atof(strtok(new_line, " "));
        natom = atof(strtok(new_line, " "));
        fion0 = atof(strtok(new_line, " "));
        nion0 = atof(strtok(new_line, " "));
        cellID = atoi(strtok(new_line, " "));
        */
        x[i] = xLoc;
        y[i] = yLoc;
        z[i] = zLoc;
        vx[i] = velx;
        vy[i] = vely;
        vz[i] = velz;
        ndencell[i] = nion0;
        temp[i] = t;
        fion[i] = fion0;
        cellnum[i] = cellID;
        zmfrac[i] = snIa + snII;
        Lcell[i] = size;
        i++;
    }

    fclose(fp);
}




double getamu(char *tranilist, char *ionlabel){

    int iflag, k, j, i;
    double dum, imass, mamu;
    char name[50], ilabel[50], tlabel[50];
    char new_line[500];
    FILE *fp = fopen(tranilist, "r");
    if (fp==NULL){
        printf("\n\nERROR in getamu in files7.c\n");
        printf("Cannot open %s\n", tranilist);
        printf("Exiting...\n");
        exit(1);
    }
    
    mamu = 0.0;

    // Loop past header
    fgets(new_line, sizeof(new_line), fp);
    
    // Loop through the file looking for the ion being used
    while(fgets(new_line, sizeof(new_line), fp)){
    
        sscanf(new_line, "%d %s %d %d %s %s %lf %lf %lf %lf", &iflag, name, &k, &j, ilabel, tlabel, &dum, &dum, &dum, &imass);
        if (strcmp(ilabel, ionlabel)!=0 && iflag==1){
            mamu = imass;
        }
    }
    
    fclose(fp);
    
    // Check that the ion was found
    if (mamu==0.0){
        printf("ERROR(getamu): no ion match in %s", tranilist);
        exit(1);
    }

    return mamu;
}





// make the file names for the lines files; similar to "cutfname"
// above, which is better commented
void mkfname(char *infile, char *galID, char *ion, char *lostag, char *linesfile){

    const char delim[2] = ".";
    char *input;
    // Format of infile:
    // galID.ion.los####.dat

        
    // Format of linesfile
    // galID.ion.los####.lines
    strcpy(linesfile, galID);
    strcat(linesfile, ".");
    strcat(linesfile, ion);
    strcat(linesfile, ".");
    strcat(linesfile, lostag);
    strcat(linesfile, ".lines");
}






// write the data to the .lines files, which will be used by
// specsynth to create the spectra
void wrtlines(double zgal, double *zline, double *Nline, double *bline, int *cellnum, int ndata){

    int i;
    
    char outfile[] = "linesout.lines";
    FILE *fp = fopen(outfile, "w");
    fprintf(fp, "%lf\n", zgal):
    for (i=0; i<ndata; i++){
        if (log10(Nlines[i]>9.0)){
            fprintf(fp, "%lf \t %lf \t %lf \t %d \n", zline[i], log10(Nline[i]), bline[i], cellnum[i]);
        }
    }

}






// write the data to the *.losdata  files
void wrtlosdata( double Slos, double Rgal, double zline, double vlos, double vabs, double dlos, double ndencell, double fion, double zmfrac, double Nline, double temp, double bline, double Vgalt, double vrp, double V_theta, double V_phi, double vzp, double xp, double yp, double zp, double rp, double theta, double phi, int cellnum){

    int i;

    FILE *fp = fopen(unitlosfile, "w");
    
    fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", Slos, Rgal, zline, vlos, log10(dlos), log10(ndencell), log10(fion), log10(amfrac), log10(Nline), log10(temp), bline, Vgalt, vrp, V_theta, V_phi, vap , xp, yp, zp, rp, theta, phi, crp/Vgalt, V_theta/Vgalt, V_phi/Vgalt, vzp/Vgalt);

    fclose(fp);
    
}





