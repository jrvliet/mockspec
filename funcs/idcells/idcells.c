
/* 
    Reads in the list of celsl found by cellfinder.c
    Finds the cells in the ion boxes
    Passes it to los7
*/

#include "idcells-subs.h"
#include "los7.h"
#include "datatypes.h"
#include "files7.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Run per ion
int idcells(char *ion, struct galaxy gal, int numcores, int losnum, 
            struct orient oprops){

    int i, numCells, cellnum, index;
    double losB, losPhi, mamu, zgal;
    char boxfile[100], new_line[1000];
    struct los losprops;
    FILE *cellsfp, *losfp, *linesfp;
    
    char losname[10], linesfile[80], losdata[80];
    char tranilist[80] = "Mockspec.transisions";
    
    // Get the mass of the ion
    mamu = getamu(tranilist, ion);
    

    // Read in the LOS properties, stored in lines.info file
    
    // Read in the gas box
    build_box_name(gal, ion, boxfile);
    
    // Get the number of cells in the box
    numCells = box_size(boxfile);

    // Allocate memory for the box file
    /*
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
    */
    double size, x, y, z, vx, vy, vz, dense, temp;
    double snII, snIa, natom, fion, nion, alpha, zmet;
    unsigned long cellid;
    double tph, trec, tcoll, tcool;

//    struct cell allcells[numCells];
    struct cell *allcells0 = (struct cell *)calloc(numCells, 
                                            sizeof(struct cell));
    struct cell *allcells = allcells0;

    FILE *fp = fopen(boxfile, "r");
    i = 0;
    while( fgets( new_line, sizeof(new_line), fp) ){
        
        sscanf(new_line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf" 
                         "%lf %lf %lf %lu %lf %lf %lf %lf", 
                          &size, &x, &y, &z, &vx, &vy, &vz, &dense, &temp, 
                          &snII, &snIa, &natom, &fion, &nion, &alpha, &zmet, 
                          &cellid, &tph, &trec, &tcoll, &tcool);

        allcells->size = size;
        allcells->x = x;
        allcells->y = y;
        allcells->z = z;
        allcells->vx = vx;
        allcells->vy = vy;
        allcells->vz = vz;
        allcells->dense = dense;
        allcells->temp = temp;
        allcells->snII = snII;
        allcells->snIa = snIa;
        allcells->natom = natom;
        allcells->fion = fion;
        allcells->nion = nion;
        allcells->alpha = alpha;
        allcells->zmet = zmet;
        allcells->cellid = cellid;
        i++;
/*        size++;        x++;        y++;        z++;        
vx++;        vy++;        vz++;        dense++;        temp++;        snII++;
snIa++;        natom++;        fion++;        nion++;        alpha++;        
zmet++;        cellid++;        tph++;        trec++;        tcoll++;tcool++;*/
    }

    fclose(fp);
    allcells = allcells0;

    // Loop over all lines of sight
    for (i=1; i<=losnum; i++){

        // Generate the .losdata file and open it
        mkfname(gal, ion, losname, linesfile);
        open_losdat(losdata, linesfile, losfp, linesfp);
        sprintf(losname, "%04d", losnum);

        // Write the first line to the .linesfile
        fprintf(linesfp, "%4.2lf\n", zgal);
        



        // Get the LOS props
        losprops = los_props(i);
    
        // Open the list of cells
        cellsfp = open_cell_list(i);
        
        // Read past the header
        fgets(new_line, sizeof(new_line), fp);
    
        // Loop over the cells
        while(fgets(new_line, sizeof(new_line), fp)){
                
            sscanf(new_line, "%d", &cellnum);
                
            // Cell numnber corresponds to the line of the gas file
            index = cellnum-1;

            // Pass into los7 to determine if it should be included
            los7(ion, losname, mamu, losprops, allcells[index], oprops, gal, linesfp);
        }
    
        // Close all files
        fclose(cellsfp);
        fclose(linesfp);
        fclose(losfp);
    
    }        

    free(allcells0);

    return 0;
}



int main(int argc, char *argv[]){


    char ion = "CIV";
    int losnum = 1;
    struct galaxy gal; 
    struct orient oprops;
    idcells(ion, gal, numcores, losnum, oprops);
    return 0;
}
