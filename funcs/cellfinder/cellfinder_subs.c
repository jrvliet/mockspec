#include "cellfinder_subs.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

void read_control_file(FILE *propfp, char *gasfile, char *galID, char *rootname,
                         double *aexpn, char *summaryLoc){

  char new_line[100], dum[100], gas_dum[100];
  int i;

  for (i=0; i<5; i++){
    fgets(new_line,sizeof(new_line),propfp);
  }

  // Get gasfile name
  sscanf(new_line, "%s %s", gasfile, dum);
  strcpy(gas_dum, gasfile);
  strtok(gas_dum, "_");
  sscanf(gas_dum, "%s", galID);
  //printf("gasfile: %s\n", gasfile);
  //printf("galID: %s\n", galID);

  // Get rootname
  fgets(new_line,sizeof(new_line),propfp);
  sscanf(new_line, "%s %s", rootname, dum);

  // Get the expansion factor
  fgets(new_line,sizeof(new_line),propfp);
  sscanf(new_line, "%lf %s", aexpn, dum);

  for (i=0; i<3; i++){
    fgets(new_line,sizeof(new_line),propfp);
  }
  
  // Get the location of the summary file
  sscanf(new_line, "%s %s", summaryLoc, dum);
  
}




void read_summary(char *galID, double *aexpn, char *cwd, double *mvir, 
                    double *rvir, double *a11, double *a12, double *a13, 
                    double *a21, double *a22, double *a23, double *a31, 
                    double *a32, double *a33, double *vpec_x, double *vpec_y, 
                    double *vpec_z){
  
  char new_line[1000];
  char a[10];
  int found;
  double expn, z;

  sprintf(a, "%.3f", *aexpn);
  // Build the filename of location
  char location[100];
  location[0] = '\0';
  strcat(location, cwd);
  strcat(location, "/../rotmat_a");
  strcat(location, a);
  strcat(location, ".txt");
  
  //printf("Summary location in read_summary: %s\n", location);

  // Open summary file
  FILE *fp = fopen(location, "r");
  if (fp==NULL){
    printf("\nERROR in read_summary function in cellfinder_subs.c\n");
    printf("Cannot open file %s \n", location);
    exit(1);
  }
  // Read past header
  fgets(new_line,sizeof(new_line),fp);

  // Loop through the file
  found = 0;
  while(fgets(new_line,sizeof(new_line),fp) && found==0) {
    sscanf(new_line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", 
            &expn, &z, mvir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33,
            vpec_x, vpec_y, vpec_z);

    if (fabs(expn-*aexpn)<=0.002){
      found = 1;
    }
  }

  if (found == 0){
    printf("Expn %0.3lf not found in summary file %s\n",*aexpn, location);
    exit(0);
  }
  
  fclose(fp);
  

}


const char *filegen(char *rootname, int losnum){
//void filegen(char *rootname, int losnum, char *outfile){
  // Generate a specific filename corresponding to each line of sight
  
    printf("\nIn filegen");
    printf("Rootname: %s\n",rootname);
    printf("losnum: %d\n",losnum);
  //char newstr[100], numstr[100];

    char numstr[100];
    char * newstr;
    newstr = malloc(sizeof(char)*100);
  
  newstr[0] = '\0';
  sprintf(numstr, "%.4d", losnum);
    //printf("Numstr: %s\n",numstr);

  strcat(newstr, "los");
    //printf("Newstr: %s\n",newstr);
  strcat(newstr, numstr);
    //printf("Newstr: %s\n",newstr);
  strcat(newstr, ".cellID.dat");
    //printf("Newstr: %s\n",newstr);

  return newstr;

}






// Write header information to output file
void write_OutfileHdr(FILE *outfp, double *aexpn, double *R0, double *phi, double *l, double *b, double *xen, double *yen, double *zen, double *losx, double *losy, double *losz, double *a11, double *a12, double *a13, double *a21, double *a22, double *a23, double *a31, double *a32, double *a33, double *Xcom, double *Ycom, double *Zcom, double *VXcom, double *VYcom, double *VZcom, double *x0, double *y0, double *z0, double *vx_obs, double *vy_obs,double *vz_obs){


  fprintf(outfp, "# Cell ID\n");

}






// Write LOS props to file
void write_LOSprops(FILE *propsfp, int losnum, double *aexpn, double *R0, double *phi, double *l, double *b, double *xen, double *yen, double *zen, double *losx, double *losy, double *losz, double *a11, double *a12, double *a13, double *a21, double *a22, double *a23, double *a31, double *a32, double *a33, double *Xcom, double *Ycom, double *Zcom, double *VXcom, double *VYcom, double *VZcom, double *x0, double *y0, double *z0, double *vx_obs, double *vy_obs,double *vz_obs){

  fprintf(propsfp, "%d \t %1.2lf \t %1.2lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \n", losnum, *R0, *phi, *xen, *yen, *zen, *losx, *losy, *losz, *a11, *a12, *a13, *a21, *a22, *a23, *a31, *a32, *a33, *Xcom, *Ycom, *Zcom, *VXcom, *VYcom, *VZcom, *x0, *y0, *z0, *vx_obs, *vy_obs, *vz_obs);
}

// Write LOS props to file
void write_LOSprops2(FILE *propsfp, double *arr){
    fprintf(propsfp, "%1.0lf \t %1.2lf \t %1.2lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.7lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \t %1.6lf \n", arr[0],  arr[1],  arr[2],  arr[3],  arr[4],  arr[5],  arr[6],  arr[7],  arr[8],  arr[9],  arr[10],  arr[11],  arr[12],  arr[13],  arr[14],  arr[15],  arr[16],  arr[17],  arr[18],  arr[19],  arr[20],  arr[21],  arr[22],  arr[23],  arr[24],  arr[25],  arr[26],  arr[27],  arr[28],  arr[29]);

}




