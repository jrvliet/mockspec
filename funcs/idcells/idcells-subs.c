

#include "idcells-subs.h"
#include "datatypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



void build_box_name( struct galaxy gal, char *ion, char *boxfile){

    // File name format:
    // <galID>._GZa<expn>.<ion>.txt
    strcpy(boxfile, gal.galID);
    strcat(boxfile, "_GZa");
    strcat(boxfile, gal.expn);
    strcat(boxfile, ".");
    strcat(boxfile, ion);
    strcat(boxfile, ".txt");
}


int box_size( char *boxfile ){

    int count = 0;
    char new_line[1000];
    FILE *fp = fopen(boxfile, "r");
    while(fgets(new_line, sizeof(new_line), fp)){
        count++;
    }
    fclose(fp);
    return count-2;  // Subtract off header
}



struct los los_props(int losnum){

    struct los props;
    char new_line[1000];
    int found = 0;
    int num;
    double b, phi, xen, yen, zen;
    double losx, losy, losz;
    double a11, a12, a13, a21, a22, a23, a31, a32, a33;
    double xcom, ycom, zcom, vxcom, vycom, vzcom;
    double x0, y0, z0, vxobs, vyobs, vzobs;
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
            props.phi = phi;
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
    fclose(fp);

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



FILE* open_cell_list(int num){

    char filename[200], losnum[10];
    
    // Filename structure
    //    los<losnum>.cellID.dat
    sprintf(losnum, "%04d", num);

    strcpy(filename, "los");
    strcat(filename, losnum);   
    strcat(filename, ".cellID.dat");
    
    FILE *fp = fopen(filename, "r");
    
    return fp;
}




void open_losdat(char *losdata, char *linesfile, FILE *losfp, FILE *linesfp){

    losfp = fopen(losdata, "w");
    linesfp = fopen(linesfile, "w");

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

}

