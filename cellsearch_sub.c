
#include "cellsearch_sub.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void cellsearch(int nlos, char *GZfile, double *x0, double *y0, double *z0, double *lb, double *mb, double *nb, double *xen, double *yen, double *zen, double *xex, double *yex, double *zex, double *lowervel, double *uppervel, FILE *fpout, FILE *fplog, double *a11, double *a12, double *a13, double *a21, double *a22, double *a23, double *a31, double *a32, double *a33){


  //  int width = 15;
  int error;
  int nhdr_cube = 2;
  int i, j, k;
  int jmax, kmax;
  int in_cell, RowN, onWall;
  char new_line[1000], str[50];
  double p1_x, p1_y, p1_z;
  double p2_x, p2_y, p2_z;
  double rx_xpt, ry_xpt, rz_xpt;
  double Dx, Dy, Dz;
  double dmin, r_max, radvel;

  // Columns in gas file
  double cellsize;
  double x, y, z, vx, vy, vz;
  double c8, c9, c10, c11, c12, c13, c14;

  p1_x = *xen;
  p1_y = *yen;
  p1_z = *zen;
  p2_x = *xex;
  p2_y = *yex;
  p2_z = *zex;
  Dx = *lb;
  Dy = *mb;
  Dz = *nb;

  // Max possible length of los vector, or distance between observer
  // and point of exit through box
  r_max = sqrt( pow(p2_x-p1_x, 2.0) + pow(p2_y-p1_y, 2.0) + pow(p2_z-p1_z, 2.0) );

  // The gas file (GZfile) used to be opened here, but I now open the 
  // gas file in the main program so that I can use the header lines
  // in the output file


  FILE *fpgas = fopen(GZfile, "r");
  //  printf("Opened GZ file: %s\n", GZfile);

  // Read past header
  //  printf("About to read past header of size %d\n",nhdr_cube);
  for (i=0; i<nhdr_cube; i++){
    fgets(new_line, sizeof(new_line), fpgas);
    //    printf("Line: %s\n",new_line);
  }

  //  printf("Read past header\n");
  // Scan through all the gas cells in the box. Have each cell check and 
  // report back if a segment of the line-of-sight vector passes
  // through their borders

  j = 0;
  k = 0;
  kmax = 0;      // Number of degereate detections in GZ file
  RowN = 0;      // Row number of cell in GZ file, which serves as cell ID
  in_cell = 0;
  onWall  = 0;

  while(fgets(new_line, sizeof(new_line), fpgas)){
    
    //    printf("In while loop in cellsearch\n");

    RowN++;

    sscanf(new_line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &cellsize, &x, &y, &z, &vx, &vy, &vz, &c8, &c9, &c10, &c11, &c12, &c13, &c14);

    // Old Rates files had coordinates in pc, new ones are in kpc
    // However, this version works off of the direct simulation output
    // and not the rates outputs, the coordinates are in pc
    // These need to be converted to kpc
    cellsize = cellsize/1000.0;

    // Get the minimum distance between the center of a gas cell
    // and the parametric line describing the LOS vector
    dmin = get_dmin(cellsize, x, y, z, p1_x, p1_y, p1_z, p2_x, p2_y, p2_z);

    //    printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", x, y, z, p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, dmin);

    // Get the cell's radial velocity
    radvel = get_radvel(x, y, z, vx, vy, vz);
    
    if (dmin<cellsize*0.5*sqrt(3.0) && radvel<=*uppervel && radvel>=*lowervel){
      //      if (dmin>20.0){
      //	printf("dmin: %lf \t cellsize: %lf \t Test: %lf \n",dmin, cellsize, cellsize*0.5*sqrt(3.0));
      //      }
      // If the shortest, perpendicular distance to the parameteric line is less
      // than half the diagonal distance through the cell (which is the maximum
      // distance you can be from the cell center and yet still be within the cell
      // borders), AND
      // the radial velocity of the cell is within the bounds of the limits THEN
      // flag the cell

      
      // Need confirmation that line does actually pass thru the cell
      confirm_cell(cellsize, x, y, z, dmin, p1_x, p1_y, p1_z, p2_x, p2_y, p2_z, Dx, Dy, Dz, r_max, *x0, *y0, *z0, &in_cell, &onWall, &rx_xpt, &ry_xpt, &rz_xpt);
      //      printf("Passed confirm cell\n");

      if (in_cell==1){
	j+=1;
	sprintf(str,"%d\n",RowN);
	//	printf("%p \t %s",fpout, str);
	error = fprintf(fpout, "%s", str);
      }

      if (onWall==1){
	k+=1;
	sprintf(str,"%d\n",RowN);
	printf("%p\n",fplog);
	fprintf(fplog, "%s\n", str);
      }

    }
  }

  jmax = j;
  kmax = k;
  
  // If k=0, then no degenerate detections were found, and all cells
  // along the LOS should be accounted for in "loscells.dat"
  /*
  if (k==0){
    printf("Number of los cells = %d\n",jmax);
  }
  else {
    printf("Number of degenerate detections = %d\n", kmax);
  }
  */
  fclose(fpgas);

}






double get_dmin(double cellsize, double x, double y, double z, double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z){

  
  // Calculate the minimum distance between the center of a gas cell 
  // and the parametric line describing the line-of-sight vector.
  // 
  // The equation for the minimum distance was taken from:
  // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
  //     
  // Let i,j,k denote unit vectors along the x,y,z axes.
  //
  // Let r be a vector from the origin to point (x,y,z):     r  = xi + yj + zk
  // Let r1 be a vector from the origin to point (p1_x,p1_y,z1): r1 = p1_x i + p1_y j + p1_z k
  // Let r2 be a vector from the origin to point (p2_x,p2_y,p2_z): r2 = p2_x i + p2_y j + p2_z k
  // The points (p1_x,p1_y,p1_z) and (p2_x,p2_y,p2_z) lie on a parametric line described by the vector equation: r0 + Dt, where 
  //   r0 is the origin,
  //   D represents the direction cosines Dx,Dy,Dz
  //   t is a parameter.
  //
  // The shortest distance between a point (x,y,z) and a line described by two points (p1_x,p1_y,p1_z) and (p2_x,p2_y,p2_z) is:
  //
  //              |(r-r1) X (r-r2)|
  //      dmin =  -----------------
  //                  |(r2-r1)|
  //
  // Since the origin is a point on the parametric line, r1 = 0i + 0j + 0k.

  double dmin = 0.0;
  double dx1, dy1, dz1, dx2, dy2, dz2;
  double term1, term2, term3;
  double mag_xprod, mag_x2;
  
  dx1 = x - p1_x;
  dy1 = y - p1_y;
  dz1 = z - p1_z;
  dx2 = x - p2_x;
  dy2 = y - p2_y;
  dz2 = z - p2_z;

  // Terms from taking the cross product of (r-r1)x(r-r2)
  term1 = (dy1*dz2 - dz1*dy2);
  term2 = (dx1*dz2 - dz1*dx2);
  term3 = (dx1*dy2 - dy1*dx2);

  // Magnitude of the cross product |(r-r1) X (r-r2)|
  mag_xprod = sqrt(term1*term1 + term2*term2 + term3*term3);
    
  // Magnitude of (r2-r1) = |(r2-r1)| = |(r2-0)| = |r2|
  mag_x2 = sqrt( pow((p2_x-p1_x),2) + pow((p2_y-p1_y),2) + pow((p2_z-p1_z),2));

  dmin = mag_xprod / mag_x2;
  
  return dmin;

}


 






void confirm_cell(double cellsize, double x, double y, double z, double dmin, double p1_x, double p1_y, double p1_z, double p2_x, double p2_y, double p2_z, double Dx, double Dy, double Dz, double r_max, double x0, double y0, double z0, int *in_cell, int *onWall, double *rx_xpt, double *ry_xpt, double *rz_xpt){

  //     Confirm that the parametric line passes thru the cell borders by calculating the (x,y,z) positions
  //     of the cell corner.  Calculate the point (x,y,z) at which the perpendicular line dmin intersects
  //     the parametric line.  See if the intersection point lies within or without the cell borders.
  //
  //     Crude schematic:
  //
  //                                  * (p2_x,p2_y,p2_z)=los vector end point
  //                                  |
  //                                  |  dmin
  //              intersection pt --> |---------* (x,y,z)=cell center
  //                                  |
  //                                  |
  //                   los vector --> |
  //                                  * (p1_x,p1_y,p1_z)=los vector starting point


  double thresh = 1.0e-6;
  *in_cell = 0;
  *onWall = 0;
  double cellHL = 0.5 * cellsize;
  double xLeft, xRight, yFront, yBack, zTop, zBottom;
  double d_x1, r_xpt, rdist;
  double r_xpt_obs;
  double xLmin, xLmax, xRmin, xRmax;
  double yFmin, yFmax, yBmin, yBmax;
  double zTmin, zTmax, zBmin, zBmax;
  

  // Equations of all 6 cell walls:
  xLeft   = x-cellHL;
  xRight  = x+cellHL;
  yFront  = y-cellHL;
  yBack   = y+cellHL;
  zBottom = z-cellHL;
  zTop    = z+cellHL;

  // Calculate intersection point:
  // This is a bit redundant, as dx1 and dx2 were calculated in get_dmin() above
  // Square of distance between cell center and origin of parametric line
  d_x1 = pow((x-p1_x),2) + pow((y-p1_y),2) + pow((z-p1_z),2);

  if (fabs(d_x1 - dmin) < thresh){
    // In case the difference is a very small number
    r_xpt = 0.0;
  }
  else {
    // Distance along parametric line from point of entry to the 
    // intersection or crossing point
    r_xpt = sqrt(d_x1 - pow(dmin,2.0));
  }


  // Find vlaue of parameter "t" corresponding to the point r_xpt on the parametric line
  rdist = sqrt( pow((x-x0),2) + pow((y-y0),2) + pow((z-z0),2));   // Length of vector from observer position to cell center
  r_xpt_obs = sqrt( pow(rdist,2) - pow(dmin,2) );   // Distance from observer position to crossing point

   
  *rx_xpt = p1_x + Dx * r_xpt;  // x-coordinate of intersection point
  *ry_xpt = p1_y + Dy * r_xpt;  // y-coordinate of intersection point
  *rz_xpt = p1_z + Dz * r_xpt;  // z-coordinate of intersection point

  
  // Does the point (rx_xpt,ry_xpt,rz_xpt) lie inside or outside the cell walls?
  if (*rx_xpt>xLeft && *rx_xpt<xRight && *ry_xpt>yFront && *ry_xpt<yBack && *rz_xpt>zBottom && *rz_xpt<zTop ){
    *in_cell = 1;
  }
  
  xLmin = xLeft - thresh;
  xLmax = xLeft + thresh;
  xRmin = xRight - thresh;
  xRmax = xRight + thresh;
  
  if (*rx_xpt>xLmin && *rx_xpt<xLmax){
    if (*ry_xpt>yFront && *ry_xpt<yBack && *rz_xpt>zBottom && *rz_xpt<zTop){
      *onWall = 1;
    }
  }
  
  if (*rx_xpt>xRmin && *rx_xpt<xRmax){
    if (*ry_xpt>yFront && *ry_xpt<yBack && *rz_xpt>zBottom && *rz_xpt<zTop){
      *onWall = 1;
    }
  }

  
  yFmin = yFront - thresh;
  yFmax = yFront + thresh;
  yBmin = yBack - thresh;
  yBmax = yBack + thresh;
  
  if (*ry_xpt>yFmin && *ry_xpt<yFmax){
    if (*rx_xpt>xLeft && *rx_xpt<xRight && *rz_xpt>zBottom && *rz_xpt<zTop){
      *onWall = 1;
    }
  }
  if (*ry_xpt>yBmin && *ry_xpt<yBmax){
    if (*rx_xpt>xLeft && *rx_xpt<xRight && *rz_xpt>zBottom && *rz_xpt<zTop){
      *onWall = 1;
    }
  }


  zTmin = zTop - thresh;
  zTmax = zTop + thresh;
  zBmin = zBottom - thresh;
  zBmax = zBottom + thresh;
  
  if (*rz_xpt>zBmin && *rz_xpt<zBmax){
    if (*rx_xpt>xLeft && *rx_xpt<xRight && *ry_xpt>yFront && *ry_xpt<yBack){
      *onWall = 1;
    }
  }
  
  if (*rz_xpt>zTmin && *rz_xpt<zTmax){
    if (*rx_xpt>xLeft && *rx_xpt<xRight && *ry_xpt>yFront && *ry_xpt<yBack){
      *onWall = 1;
    }
  }

}





// Function to calcuate the radial velocity of the gas in a cell
double get_radvel(double x, double y, double z, double vx, double vy, double vz){
    
  // Turns out you don't apply the rotation matrix
  // Since we are only dealing with spherical radial velocity,
  // the rotation is unnessecary
  double lenrg, dotprod, vrad;

  lenrg = sqrt( x*x + y*y + z*z );
  dotprod = vx*x + vy*y + vz*z;
  vrad = dotprod / lenrg;

  return vrad;

}


