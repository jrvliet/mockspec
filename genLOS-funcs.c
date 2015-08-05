

int inBox(double x, double y, double z, double boxsize){

    if (x<boxsize && x>-1*boxsize && y<boxsize && y>-1*boxsize && z<boxsize && z>-1*boxsize){
        result = 1;
    }       
    else{
        result = 0;
    }
    return result;
}


void findEnds(double px, double py, double pz, double dx, double dy, double dz, double boxsize, double *xen, double *yen, double *zen, double *xex, double *yex, double *zex){

    double size = boxsize/2.0;
    double arrayStart = -5000.;
    double arrayEnd = 5000.;
    double arrayNum = 5.0e5;
    double step = (arrayEnd-arrayStart) / arrayNum;

    int i, point1In, point2In;
    int 
    double x1, y1, z1, x2, y2, z2;
    double *t;   // Steps along LOS
    t = (double *)calloc(arrayNum, sizeof(double));
    for (i=0; i<arrayNum; i++){
        t[i] = arrayStart + i*step;
    }
    
    // Loops through array to find the ends
    for (i=0; i<arrayNum-1; i++){
        x1 = px + dx*t[i];
        y1 = py + dy*t[i];
        z1 = pz + dz*t[i];
        
        x2 = px + dx*t[i+1];
        y2 = py + dy*t[i+1];
        z2 = pz + dz*t[i+1];
        
        point1In = inBox(x1, y1, z1, size);
        point2In = inBox(x2, y2, z2, size);
        
        if (point1In==0 && point2In==1){
            &xen = (x1+x2)/2.0;
            &yen = (y1+y2)/2.0;
            &zen = (z1+z2)/2.0;
            &ten = (t[i]+t[i+1])/2.0;
        }
        if (point1In==1 && point2In==2){
            &xex = (x1+x2)/2.0;
            &yex = (y1+y2)/2.0;
            &zex = (z1+z2)/2.0;
            &tex = (t[i]+t[i+1])/2.0;
        }
    }
    free(t);
}

////////////////////////////////////////////////////////////////

void read_summary(char *galID, double *aexpn, char *summaryLoc, double *mvir, double *rvir, double *a
11, double *a12, double *a13, double *a21, double *a22, double *a23, double *a31, double *a32, double
 *a33, double *vpec_x, double *vpec_y, double *vpec_z){
  char new_line[1000];
  int found;
  double expn, z;

  // Build the filename of location
  char location[100];
  location[0] = '\0';
  strcat(location, summaryLoc);
  strcat(location, galID);
  strcat(location, ".dat");

  // Open summary file
  FILE *fp = fopen(location, "r");

  printf("%s\n",location);
  // Read past header
  fgets(new_line,sizeof(new_line),fp);
  fgets(new_line,sizeof(new_line),fp);

  // Loop through the file
  found = 0;
  while(fgets(new_line,sizeof(new_line),fp) && found==0) {
    //    printf("Line in summary file\n %s\n",new_line);
    sscanf(new_line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf ", &expn, &z, mv
ir, rvir, a11, a12, a13, a21, a22, a23, a31, a32, a33, vpec_x, vpec_y, vpec_z);


    if (expn==*aexpn){
      found = 1;
    }
  }

  if (found == 0){
    printf("Expn %0.3lf not found in summary file %s\n",*aexpn, location);
    exit(0);
  }
  
  fclose(fp);
  
}







