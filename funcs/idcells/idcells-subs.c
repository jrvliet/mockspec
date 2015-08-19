

#include "idcells-subs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct los{
    int num;
    double b;
    double phi;
    double xen;
    double yen;
    double zen;
    double losx;
    double losy;
    double losz;
    double a11;
    double a12;
    double a13;
    double a21;
    double a22;
    double a23;
    double a31;
    double a32;
    double a31;
    double xcom;
    double ycom;
    double zcom;
    double vxcom;
    double vycom;
    double vzcom;
    double x0;
    double y0;
    double z0;
    double vxobs;
    double vyobs;
    double vzobs;


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

        int count = 0;c
    char new_line[1000];
    FILE *fp = fopen(boxfile, "r");
    while(fgets(new_line, sizeof(new_line), fp)){
        count++;
    }
    fcloes(fp);
    return count-2;  // Subtract off header
}



struct los losProps(int losnum){

    struct los props;
    char new_line[1000];
    int found = 0;

    FILE *fp = fopen("lines.props", "r");
    fgets(new_line, sizeof(new_line), fp);

    while(fgets(new_line, sizeof(new_line), fp)){
    
        sscanf(new_line, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
                         "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
                         "%lf %lf %lf %lf",
                          &num, &b, &phi, &xen, &yen, &zen, &losx, &losy, &losz,
                          &a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33, 
                          &xcom, &ycom, &zcom, &vxcom, &vycom, &vzcom, 
                          &x0, &y0, &z0, &vxobs, &vyobs, &vzobs);
        if (num==losnum){
            props.num = num;
            props.b = b;
            props.phi = phi
            props.xen = xen;
            props.yen = yen;
            props.zen = zen;
            props.losx = losx;
            props.losy = losy;
            props.losz = losz;
            props.a11 = a11;
            props.a12 = a12;
            props.a13 = a13;
            props.a21 = a21;
            props.a22 = a22;
            props.a23 = a23;
            props.a31 = a31;
            props.a32 = a32;
            props.a33 = a33;
            props.xcom = xcom;
            props.ycom = ycom;
            props.zcom = zcom;
            props.vxcom = vxcom;
            props.vycom = vycom;
            props.vzcom = vzcom;
            props.x0 = x0;
            props.y0 = y0;
            props.z0 = z0;
            props.vxobs = vxobs;
            props.vyobs = vyobs;
            props.vzobs = vzobs;
            found = 1;
        }
    }
    fp.close();

    if (found==0){
        printf("ERROR in idcells, function losProps\n");
        printf("Could not find losnum %d\n", losnum);
        printf("Exiting...\n");
        exit(1);
    }
    else{
        return props;
    }
}










