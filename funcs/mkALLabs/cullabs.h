
      integer           maxlos,maxcols,mxions,
     &                  mxtrans,mxreg,nregcols,Mmax
      parameter         (maxlos   = 10000,  ! max sightines in survey
     &                   maxcols  =    26,  ! max columns in losdata files
     &                   mxions   =    30,  ! max ions in survey
     &                   mxtrans  =   255,  ! max transitions in tranilist
     &                   mxreg    =    20,  ! max number of abs regions
     &                   nregcols =    15,  ! # of cols in input regabs files
     &                   Mmax     = maxlos*mxions) 
