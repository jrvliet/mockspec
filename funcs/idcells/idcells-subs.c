

#include "idcells-subs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void buildBoxName( char *galID, char *expn, char *ion, char *boxfile){

    // File name format:
    // <galID>._GZa<expn>.<ion>.txt
    strcpy(boxfile, galID);
    strcat(boxfile, "_GZa");
    strcat(boxfile, expn);
    strcat(boxfile, ".");
    strcat(boxfile, ion);
    strcat(boxfile, ".txt");
}


int boxSize( char *boxfile ){

    int count = 0;
    char new_line[1000];
    FILE *fp = fopen(boxfile, "r");
    while(fgets(new_line, sizeof(new_line), fp)){
        count++;
    }
    fcloes(fp);
    return count-2;  // Subtract off header
}
    
