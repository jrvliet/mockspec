//!/usr/bin/python
//
// Filename: cellfinder.py
//
// Purpose: 
//  This is a version of FBcellfinder but does not 
//  This does not get any information from the name of the galaxy file
//  All information is pulled from the galaxy props file
//  This is to allow the use of galaxies with naming conventions other than
//  Ceverino's. 
//
//  Run with:
//    python FBcellfinder_v7.py <Gal Props File> <lower bound> <upper bound>
//  
//  Example:
//    python FBcellfinder_v7.py gal_props.dat 0 inf
//

#include "cellfinder-subs.h"
#include "lmn_cellfinder.h"
#include "cellsearch_sub.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



int main(int argc, char *argv[]){



  FILE *propfp0 = fopen("gal_props.dat", "r");
  FILE *propfp = propfp0;
  double lowvel, upvel;

  upvel = 1.0e30;
  lowvel = -1.9e30;
  

  int nlos, losnum;
  char new_line[100], new_line2[100];
  char gasfile[100], galID[100], rootname[100], summaryLoc[100];
  char outfile[100];
  FILE *outfp0, *outfp;
  double aexpn, mvir, rvir;
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;
  double vpec_x, vpec_y, vpec_z;
  double Xcom, Ycom, Zcom;
  double VXcom, VYcom, VZcom;
  double l, b;
  double vx_obs, vy_obs, vz_obs;
  double boxlength;
  double losx, losy, losz, x0, y0, z0;  // Describe unit vector
  double xlen, ylen, zlen, loslength;

  // Read input parameter file which contains all of the galaxy information
  read_control_file(propfp, gasfile, galID, rootname, &aexpn, summaryLoc);

  // Read summary file
  read_summary(galID, &aexpn, summaryLoc, &mvir, &rvir, &a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33, &vpec_x, &vpec_y, &vpec_z);

  Xcom = 0.0;
  Ycom = 0.0;
  Zcom = 0.0;
  VXcom = 0.0;
  VYcom = 0.0;
  VZcom = 0.0;
  l = 0.0;
  b = 90.0;
  vx_obs = 0.0;
  vy_obs = 0.0;
  vz_obs = 0.0;
  
  // Write the ouput log file header
  FILE *logfp = fopen("CellWalls.log", "w");
  fprintf(logfp, "The LOS vector runs along the cell walls of the following cells:\n");
  fprintf(logfp, "LOS # \t Cell ID #\n");

  // Get information from gasfilename
  boxlength = 4.0*rvir;

  // Start looping over lines of sight
  // Print header to screenfp
  FILE *screenfp = fopen("cellfinder.screen", "w"); 
  fprintf(screenfp, "===============================================\n");
  fprintf(screenfp, "Gas file: %s\n",gasfile);
  fprintf(screenfp, "===============================================\n");

  // Read in los of sight from posnlist.dat
  FILE *losfp = fopen("lines.dat", "r");
  FILE *infofp = fopen("lines.info", "r");
  FILE *propsfp = fopen("lines.props", "w");

  // Get number of lines of sight
  nlos = 0;
  while(fgets(new_line,sizeof(new_line),losfp)){
    if(new_line[0] != '#')
      nlos++;
  }
  // Point fpin back to the beginning of the file
  fseek(losfp,0L,SEEK_SET);

  fprintf(screenfp, "Number of lines of sight: %d\n",nlos);

  // Write header to lines.props file
  fprintf(propsfp, "LOS \t b \t %-6s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \t %-10s \n", "phi", "xen", "yen", "zen", "losx", "losy", "losz", "a11", "a12", "a13", "a21", "a22", "a23", "a31", "a32", "a33", "Xcom", "Ycom", "Zcom", "VXcom", "VYcom", "VZcom", "x0", "y0", "z0", "vx_obs", "vy_obs", "vz_obs");


  double xen, yen, zen, xex, yex, zex;
  double R0, phi, incline;

  losnum = 0;
  while(fgets(new_line,sizeof(new_line),losfp)){
    fgets(new_line2,sizeof(new_line),infofp);
    if(new_line[0] != '#'){
      
      // Read in line from lines.dat
      sscanf(new_line, "%lf %lf %lf %lf %lf %lf", &xen, &yen, &zen, &xex, &yex, &zex); 
      // Read in line from lines.info
      sscanf(new_line2, "%d %lf %lf %lf", &losnum, &R0, &phi, &incline); 
      // Get the direction cosines (losx, losy, losz) and unit vector
      // describing los vector through box (tail at entrance point, headt at (x0, y0, z0)
      lmn4(xen, yen, zen, xex, yex, zex, Xcom, Ycom, Zcom, &losx, &losy, &losz, &x0, &y0, &z0);
      fprintf(screenfp, "\nPOEntry: %3.6lf %3.6lf %3.6lf\n",xen,yen,zen);
      fprintf(screenfp, "POExit:   %3.6lf %3.6lf %3.6lf\n",xex,yex,zex);
      
      xlen = xex-xen;
      ylen = yex-yen;
      zlen = zex-zen;
      loslength = sqrt(xlen*xlen + ylen*ylen + zlen*zlen);
      fprintf(screenfp, "LOS Length: %3.2lf\n",loslength);

      // Generate a unique filename corresponding to current los in loslist
      strcpy(outfile, filegen(rootname, losnum));
      fprintf(screenfp, "Output file: %s\n",outfile);
      
      // Open outfile
      outfp0 = fopen(outfile, "w");
      outfp = outfp0;

      // Write header to output file
      write_OutfileHdr(outfp, &aexpn, &R0, &phi, &l, &b, &xen, &yen, &zen, &losx, &losy, &losz, &a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33, &Xcom, &Ycom, &Zcom, &VXcom, &VYcom, &VZcom, &x0, &y0, &z0, &vx_obs, &vy_obs, &vz_obs);

      // Write LOS properties to file
      write_LOSprops(propsfp, losnum, &aexpn, &R0, &phi, &l, &b, &xen, &yen, &zen, &losx, &losy, &losz, &a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33, &Xcom, &Ycom, &Zcom, &VXcom, &VYcom, &VZcom, &x0, &y0, &z0, &vx_obs, &vy_obs, &vz_obs);

      outfp = outfp0;
      // Find all the cells that lie along the line of sight as 
      // defined by the points of entry and exit into and out of the box
      cellsearch(nlos, gasfile, &x0, &y0, &z0, &losx, &losy, &losz, &xen, &yen, &zen, &xex, &yex, &zex, &lowvel, &upvel, outfp, logfp, &a11, &a12, &a13, &a21, &a22, &a23, &a31, &a32, &a33);
      fclose(outfp0);
    }

  }

  fclose(propsfp);
  fclose(screenfp);
  return 0;
}
