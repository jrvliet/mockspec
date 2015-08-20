
/* 
    Reads in the list of celsl found by cellfinder.c
    Finds the cells in the ion boxes
    Passes it to los7
*/

#include "idcells-subs.h"
#include "los7.h"
#include <stdio.h>1
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Run per ion
int idcells(char *ion, char *galID, char *expn, int numcores, int losnum){

    int i, numCells;
    double losB, losPhi, mamau;
    char boxfile[100], new_line[1000];
    struct los props;
    FILE *cellsfp;

    char tranilist[80] = "Mockspec.transisions";

    // Get the mass of the ion
    mamu = getamu(tranilist, ion);
    

    // Read in the LOS properties, stored in lines.info file
    
    // Read in the gas box
    buildBoxName(galID, expn, ion, boxfile);
    
    // Get the number of cells in the box
    numCells = boxSize(boxfile);

    // Allocate memory for the box file
    double *size = (double *)calloc(numCells, sizeof(double));
    double *x = (double *)calloc(numCells, sizeof(double));
    double *y = (double *)calloc(numCells, sizeof(double));
    double *z = (double *)calloc(numCells, sizeof(double));
    double *vx = (double *)calloc(numCells, sizeof(double));
    double *vy = (double *)calloc(numCells, sizeof(double));
    double *vz = (double *)calloc(numCells, sizeof(double));
    double *dense = (double *)calloc(numCells, sizeof(double));
    double *temp = (double *)calloc(numCells, sizeof(double));
    double *snII = (double *)calloc(numCells, sizeof(double));
    double *snIa = (double *)calloc(numCells, sizeof(double));
    double *natom = (double *)calloc(numCells, sizeof(double));
    double *fion = (double *)calloc(numCells, sizeof(double));
    double *nion = (double *)calloc(numCells, sizeof(double));
    double *alpha = (double *)calloc(numCells, sizeof(double));
    double *zmet = (double *)calloc(numCells, sizeof(double));
    unsigned long *cellid = (unsigned long *)calloc(numCells, 
                                             sizeof(unsigned long));
    double *tph = (double *)calloc(numCells, sizeof(double));
    double *trec = (double *)calloc(numCells, sizeof(double));
    double *tcoll = (double *)calloc(numCells, sizeof(double));
    double *tcool = (double *)calloc(numCells, sizeof(double));
    
    FILE *fp = fopen(boxfile, "r");

    while( fgets( new_line, sizeof(new_line), fp) ){
        
        sscanf(new_line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf" 
                         "%lf %lf %lf %lu %lf %lf %lf %lf", 
                          size, x, y, z, vx, vy, vz, dense, temp, snII, snIa,
                          natom, fion, nion, alpha, zmet, cellid, 
                          tph, trec, tcoll, tcool);
        size++;
        x++;
        y++;
        z++;
        vx++;
        vy++;
        vz++;
        dense++;
        temp++;
        snII++;
        snIa++;
        natom++;
        fion++;
        nion++;
        alpha++;
        zmet++;
        cellid++;
        tph++;
        trec++;
        tcoll++;
        tcool++;
    }

    fp.close()


    // Loop over all lines of sight
    for (i=1; i<=losnum; i++){

        // Generate the .losdata file and open it
        mkfname(losdata, galID< ion, losname, linesfile);
        open_losdat(losdata, linesfile, losfp, linesfp);



        // Get the LOS props
        props = losProps(i);
    
        // Open the list of cells
        cellsfp = open_cell_list(i);
        
        // Read past the header
        fgets(new_line, sizeof(new_line), fp);
    
        // Loop over the cells
        while(fgets(new_line, sizeof(new_line), fp){
                
            sscanf(new_line, "%d", &cellnum);
                
            // Cell numnber corresponds to the line of the gas file
            index = cellnum-1;

            // Pass into los7 to determine if it should be included


    return 0;
}

