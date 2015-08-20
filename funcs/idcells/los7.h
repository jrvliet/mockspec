#ifndef LOS7_H
#define LOS7_H

#include <stdio.h>
#include "datatypes.h"


int los7(char *ion, char *losname, double mamu, struct los losprops, 
            struct cell cprops, struct orient oprops, struct galaxy gprops,
            FILE *linesfp);


#endif
