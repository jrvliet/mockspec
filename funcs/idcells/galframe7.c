
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


//     the first routine performs matrix rotation to obtain the cell
//     cartesian positions into the reference frame of the galaxy and the
//     second routine computes the spherical coordinate positions and the
//     spherical coordinate components of the cell velocity in the galaxy
//     reference frame


void rotate(double a[][3], double x, double y, double z, double *xp, double *yp, double *zp){

    // rotate matrix operation, quantities sent down (x,y,z) have already
    // been translated to galaxy center; returned quantities are the
    // position (xp,yp,zp) in the galaxy reference frame

    *xp = a[0][0]*x + a[0][1]*y + a[0][2]*z;
    *yp = a[1][0]*x + a[1][1]*y + a[1][2]*z;
    *zp = a[2][0]*x + a[2][1]*y + a[2][2]*z;

}




void sphvels( double xp, double yp, double zp, double *rp, double *theta, double *phi, double vxp, double vyp, double vzp, double *vrp, double *vtheta, double *vphi){

    // spherical coordinate velocities, quantities sent down (xp,yp,zp)
    // were determined in routine rotate

    // theta = polar  angle    - (great circle component, streams, etc.)
    // phi   = azimuthal angle - (angular momentum component, rotation, etc)

    double ihat, jhat, khat;
    double rho;
    double vth_top, vth_bot;
    double rp0, theta0, phi0, vrp0, vtheta0, vphi0;

    // Compute the cartesian unit vectors
    rp0 = sqrt( xp*xp + yp*yp + zp*zp );
    ihat = xp/rp0;
    jhat = yp/rp0;
    khat = zp/rp0;

    // The polar coordinates; the cosine of the khat unit vector
    theta0 = acos( khat );  // returned in radians
    theta0 = theta0*57.296;        // convert to degrees

    // The azimuthal coordinates
    phi0 = atan2(yp,xp);   // preserves the quadrant (-pi<=phi<=pi)
    phi0 = phi0*57.296;

    // The radial velocity
    // Cartesian velocity components dotted into the unit vectors
    vrp0 = vxp*ihat + vyp*jhat + vzp*khat;
    
    // The polar unit vector velocity 
    // If smack on the galaxy center, this is undefined and we null it
    if (rp0!=0.0){
        vth_top = zp*(vrp0/rp0) - vzp;
        vth_bot = sqrt(1.0 - pow( (zp/rp0), 2.0 ));
        vtheta0 = vth_top/vth_bot;
    }
    else{
        vtheta0 = 0.0;
    }
    
    // the azimuthal unit vector velocity; first we compute the
    // projection of the radial distance (rho) onto the galaxy plane; if
    // there is no projection (theta=+/-90) then this is undefined and we
    // null it
    rho = sqrt( xp*xp + yp*yp);
    
    if (rho!=0.0){
        vphi0 = (xp/rho) * (vyp - yp*(vxp/xp));
    }
    else{
        vphi0 = 0.0;
    }

    *rp = rp0;
    *theta = theta0;
    *phi = phi0;
    *vrp = vrp0;
    *vtheta = vtheta0;
    *vphi = vphi0;

}

    





