c for best display, widen screen until the below line does not wrap
c---------------------------------------------------------------------------
c
c     PACKAGE: Mockspec
c     PROGRAM: sysanal
c     MODULE:  main routine
c
c     DESCRIPTION
c     this code performs analysis on absorption systems for which each
c     ion/transition is stored in a ASCII file over some fixed velocity
c     range.  it is designed specifically to process output files from
c     spectrum synthesis program Mockspec/mkspectra/specsynth as part of
c     the Mockspec package
c
c
c     INCLUDED MODULES
c     access.f         homegrown for f95 (remove if f77)
c     files.f          manages file I/O
c     subs.f           various subroutines
c     features.f       routines for feature finding/manipulation
c     calc.f           routines for calculating EWs, etc.
c     spectra.f        routines for computing EW, AOD spectra, etc.
c     zabs.f           computes the systetmic redshift
c     recipes.f        various required Numerical Recipes (modified)
c     setup.f          special routine for Mock spectra analysis
c
c
c     author: Chris Churchill
c     email:  cwc@nmsu.edu
c
c     last update: May 2011
c
c
c.........................................................................
c

      PROGRAM sysanal

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      logical            detect,error
      integer            ntran(mxions)
      integer            i,j,jknt
      double precision   EWtest
      character*80       qsolist,paramlist,losnum,commflag
      character*80       tranilist,instrlist,losindex
      character*80       specname(mxions)
      character*80       trani(mxions,mxtrans),
     &                   transpecfile(mxions,mxtrans)

      include            'sysanal.com'


      error = .false.

c     grab the command line arguments

      CALL getarg(1,qsolist)     ! list of sightlines (QS0.XXX.dat files)

      IF (qsolist.eq.'') THEN        
       WRITE(6,*) 'Usage: sysanal $1 [$2] [$3] [$4] [$5]'        
       WRITE(6,*) '  $1  = filename of list of QS0.XXX.dat files'        
       WRITE(6,*) ' [$2] = survey configuration file'        
       WRITE(6,*) ' [$3] = atomic/transition data file'        
       WRITE(6,*) ' [$4] = instrument configuration file'        
       WRITE(6,*) ' '        
       WRITE(6,*) ' if $2 is null "Mockspec.runpars" is assumed.'        
       WRITE(6,*) ' it is assumed this file as 1 header line.'        
       WRITE(6,*) ' '        
       WRITE(6,*) ' if $3 is null "Mockspec.transitions" is assumed.'       
       WRITE(6,*) ' it is assumed this file has 1 header line.'        
       WRITE(6,*) ' '        
       WRITE(6,*) ' if $4 is null "Mockspec.instruments" is assumed.'       
       WRITE(6,*) ' it is assumed this file has 1 header line.'        
       WRITE(6,*) ' '        
       WRITE(6,*) ' if $5=0 ornull I/O is suppressed. All other values' 
       WRITE(6,*) ' result in I/O to screen and runlog file.'        
       WRITE(6,*) ' '        
       WRITE(6,*) ' these files need to be in the present working'        
       WRITE(6,*) ' directory: for required file formating, see'        
       WRITE(6,*) ' the help file "Mockspec/Mockspec.help"'        
       WRITE(6,*) ' '        
       STOP      
      END IF


      CALL getarg(2,paramlist)         
      If (paramlist.eq.'') paramlist = 'Mockspec.runpars'      

      CALL getarg(3,tranilist)         
      If (tranilist.eq.'') tranilist = 'Mockspec.transitions'      

      CALL getarg(4,instrlist)         
      If (instrlist.eq.'') instrlist = 'Mockspec.instruments'      

      CALL getarg(5,commflag)         
      IF ((commflag.eq.'').OR.(commflag.eq.'0')) then
        iprint = .false.
      ELSE
        iprint = .true.
      ENDIF

c     check the command line input file sto make sure they exist; if not
c     communicate and terminate

      CALL ckfiles(qsolist,paramlist,tranilist,instrlist)

c     clean the slate

      CALL clean

c
c     START
c

c     open the runlog file; communicate the input lists; configure for
c     the run parameters; read in the galaxy systemic redshift

      IF (iprint) then
       OPEN(unit=1,file='Mockspec.runlog.anal',status='unknown')
       CALL comm1(qsolist,paramlist,tranilist,instrlist)
      ENDIF

c     obtain the survey configuration

      CALL getparams(paramlist)

c     open the qsolist file (the file containing the sightlines)

      OPEN(unit=33,file=qsolist,status='old')

c     write the header for the runlog file and to screen

      IF (iprint) then
       WRITE(6,600)
       WRITE(1,600)
      ENDIF

c     loop over the qsolist 

      DO 11 i=1,nmx

c     read in the sightline filename *.dat (losnum) and write to the
c     logfile and screen

       READ(33,*,END=12) losnum
       CALL substr2(losnum,losindex)

c     jknt=1 set in routine gettransitions (left over due to historical
c     architecture of general code sysanal)

       jknt = 0
       CALL gettransitions(jknt,tranilist,losnum,specname,
     &            ntran,trani,transpecfile,losindex,error)
       IF (error) GOTO 11

c     set up the arrays and file names for the analysis

       CALL setup(jknt,ntran,trani,transpecfile,instrlist)

c     zero out all of the working arrays

       CALL zeroall

c     read in the spectra, if no cells error=.true.
c     move on to the next sightline

       error = .false.
       CALL getdata(specname(jknt),error)
       IF (error) GOTO 11

c     automated computation of the initial system master regions and
c     subregions; get the initial redshift for the auto regions

       CALL initregs

c     if no detections, the communicate no absorption found
c     move on to the next sightline

       IF (nlines.eq.0) THEN
        IF (iprint) CALL comm0(jknt,losindex,specname(jknt))
        GOTO 11
       END IF

c     compute the absorption properties

       CALL finalcalc

c     if the rest-frame EW is .ge. EWcut, then include system in the
c     survey, if not communicate rejection from the survey, account for
c     1-sigma measurement uncertainty in the cut off

C     MODIFIED to print out all LOS by commenting out below condition

       detect = .false.
       EWtest = ewtot(1)/(1.0+zbar)
       IF (EWtest.ge.EWcut(jknt)) detect= .true.
C        detect = .true.
        CALL writefiles(losnum)
C       END IF
       IF (iprint) CALL comm2(jknt,losindex,specname(jknt),detect)

 11   CONTINUE

 12   CLOSE(unit=33)
      CLOSE(unit=1)

      STOP 

 600  FORMAT(1x,3x,'LOS#',3x,'ion',8x,'Instr',7x,'EWcut',5x,'Nreg',
     &       3x,'EWr',6x,'Nsig',2x,'Note',16x,'keytran',3x,'+tran(s)')

      END

c
c-----------------------------------------------------------------------
      include 'access.f'     ! homegrown for f95 (remove if f77)
      include 'files.f'      ! manages file I/O
      include 'subs.f'       ! various subroutines
      include 'features.f'   ! routines for feature finding/manipulation
      include 'calc.f'       ! routines for calculating EWs, etc.
      include 'spectra.f'    ! routines for computing EW, AOD spectra, etc.
      include 'zabs.f'       ! computes the systetmic redshift 
      include 'recipes.f'    ! various required Numerical Recipes (modified)
      include 'setup.f'      ! special routine for Mock spectra analysis
c-----------------------------------------------------------------------
