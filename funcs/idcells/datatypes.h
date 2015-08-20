#ifndef DATATYPES_H
#define DATATYPES_H

#include <stdio.h>


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
    double a11;
    double a12;
    double a13;
    double a21;
    double a22;
    double a23;
    double a31;
    double a32;
    double a33;
};

struct orient{
    double incline;
    double losx;
    double losy;
    double losz;
};


struct cell{
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double size;
    double dense;
    double temp;
    double natom;
    double zmet;
    double alpha;
    double nion;
    double fion;
    double snII;
    double snIa;
    unsigned long cellid;
};

struct galaxy{
    char *galID;
    char *expn;
    double mvir;
    double rvir;
    double a11;
    double a12;
    double a13;
    double a21;
    double a22;
    double a23;
    double a31;
    double a32;
    double a33;
    double vxg;
    double vyg;
    double vzg;
};  


#endif
