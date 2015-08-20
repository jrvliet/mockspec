#ifndef IDCELLS_SUBS_H
#define IDCELLS_SUBS_H

#include <stdio.h>
#include "datatypes.h"

void build_box_name( struct galaxy gal, char *ion, char *boxfile);

int box_size( char *boxfile );

struct los los_props(int losnum);

FILE* open_cell_list(int num);

void open_losdat(char *losdata, char *linesfile, FILE *losfp, FILE *linesfp);

#endif
