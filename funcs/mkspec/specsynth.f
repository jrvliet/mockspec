c     PACKAGE: Mockspec
c     PROGRAM: specsynth
c     MODULE:  main program
c
c
c     DESCRIPTION 
c     this program computes the mock spectra from the output files from
c     the LOS program; it is a hacked version of the original specsynth
c     program that handles multiple transition in a single spectrum;
c     this version handles only one type of transition in a synthesized
c     spectrum (no blends)
c
c
c     INCLUDED modules:
c     const.dek
c     access.f
c     convolve.f
c     fft.f
c     instrument.f
c     voigt.f
c     model.f
c     noise.f
c     files.f
c     strings.f
c     splinecode.f
c
c
c     author: Chris Churchill
c     email:  cwc@nmsu.edu
c     
c     last update:    Sep 2015 - modified do_comm logical
c     created update: May 2011
c
c
c
c..............................................................................
c

      PROGRAM    spectra

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h' 
      include           'const.dek'
 
      logical            error
      integer            maxqsos,ntran,idum
      parameter          (maxqsos=10000)
      integer            i,j,m
      double precision   wcen
      character*80       qsolist,paramlist,losnum
      character*80       tranilist,instrlist,logflag  
      character*80       linelist(maxtrans),
     &                   trani(maxtrans)



      idum = -35
      WRITE(6,*) idum
      CALL SYSTEM_CLOCK(idum)
      WRITE(6,*) idum
      idum = -1*idum
      WRITE(6,*) idum

c     read in the list of ion names, input the ion information and book
c     keeping, input the atomic constants, etc., grab the parameters


      CALL getarg(1,qsolist)     ! list of sightlines (QS0.XXX.dat files)

      IF (qsolist.eq.'') THEN
        WRITE(6,*) 'Usage: specsynth $1 [$2] [$3] [$4] [$5]'
        WRITE(6,*) '  $1  = filename of list of QS0.XXX.dat files'
        WRITE(6,*) ' [$2] = survey configuration file'
        WRITE(6,*) ' [$3] = atomic/transition data file'
        WRITE(6,*) ' [$4] = instrument configuration file'
        WRITE(6,*) ' [$5] = logfile/screen write suppression flag'
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
        WRITE(6,*) ' if $5=0 or null, then logfile/screen communication'
        WRITE(6,*) ' is turned OFF; any other value turns it ON'
        WRITE(6,*) ' '
        WRITE(6,*) ' these files need to be in the present working'
        WRITE(6,*) ' directory: for required file formating, see'
        WRITE(6,*) ' the help file "Mockspec/Mockspec.help"'
        WRITE(6,*) ' '
        STOP
      END IF

      CALL getarg(2,paramlist)   
      IF (paramlist.eq.'') paramlist = 'Mockspec.runpars'

      CALL getarg(3,tranilist)   
      If (tranilist.eq.'') tranilist = 'Mockspec.transitions'

      CALL getarg(4,instrlist)   
      IF (instrlist.eq.'') instrlist = 'Mockspec.instruments'

      CALL getarg(5,logflag)   
      IF ((logflag.eq.'').OR.(logflag.eq.'0')) then
       iprint = .false.
      ELSE
       iprint = .true.
      ENDIF

      CALL ckfiles(qsolist,paramlist,tranilist,instrlist)

c
c     START
c

c     clean

      CALL clean


c     open the runlog file; communicate the input lists; configure for
c     the run parameters; read in the galaxy systemic redshift

      IF (iprint) then
       OPEN(unit=1,file='Mockspec.runlog.spec',status='unknown')
       CALL comm1(qsolist,paramlist,tranilist,instrlist)
      ENDIF

      CALL getparams(paramlist)  

c     open the QSOLIST file

      OPEN(unit=33,file=qsolist,status='unknown')

c
c     MAIN OUTER LOOP: loop over the qsolist 
c

      DO 11 i=1,maxqsos

c     read in the sightline filename QS0.XXX.dat (called losnum)

        READ(33,*,END=12) losnum

c     intiliza some counters 

        lines = 0
        ntran = 0

c     for this losnum file locate all the transition that we will
c     generate spectra for (see notes in routine gettransitions); on
c     return we have the number of transitions to model (ntran), the
c     character array trani(1:ntran) containing the transition names,
c     and the character array linelist(1:ntran) containing the names of
c     the files that have the reshifts, column densities, and Doppler b
c     parameters required to generate the spectra; we also have the
c     atomic constants for each transition (there are three and they are
c     global variables)

        CALL gettransitions(tranilist,losnum,ntran,trani,linelist)

c     INNER LOOP: loop over all the transition and make the spectra one
c     by one

        DO 09 j=1,ntran       

c     open the linelist file; if it DNE then skip (shouldn't ever
c     happen, it was accessed on disk in routine gettransitions); if
c     there are no lines listed (again, this should never happen) then
c     skip, communicate

          CALL readlines(linelist(j),trani(j))   
          IF (lines.eq.0) THEN
            If (iprint) CALL comm0lines(losnum,j,trani(j),linelist(j))
            GOTO 09
          ENDIF

c     assigns atomic constants to all lines for this transition and
c     obtain this transitions instrumental configuration

          CALL setatomic(j)             
          CALL configspec(instrlist,j,wcen)

c     delete any lines that are out of this spectral range; initializes
c     the spectrum continuum

          CALL droplines                 
          CALL initspectrum              

          m    = ndata

          CALL instrument(m,wcen,0)      ! the ISF FWHM over the spectrum
          CALL doabslines                ! line flux radiative transfer
          CALL dolymanlimit              ! Lyman limit break
          IF (convolving) CALL convolve  ! FFT convolution with ISF
          CALL addnoise(idum)            ! add gaussian deviates with SNR

c     write the output spectrum; communicate this transition's data and
c     the spectral configuration

          CALL output(losnum,trani(j))                

          IF (iprint) CALL comm2(losnum,j,trani(j),linelist(j))
                                             
c     go to the next transition

 09     CONTINUE

c     go to the next sightline in the qsolist

 11   CONTINUE

c     close up shop 

 12   CLOSE(unit=33)

c     write a final barline and close the runlog file

      IF (iprint) then 
       WRITE(1,*) ' ' 
       WRITE(1,*) '*******************************************'
       WRITE(6,*) ' ' 
       WRITE(6,*) '*******************************************'
       CLOSE(unit=1)
      ENDIF

      STOP
      END


      include 'access.f'
      include 'convolve.f'
      include 'fft.f'
      include 'instrument.f'
      include 'voigt.f'
      include 'model.f'
      include 'noise.f'
      include 'files.f'
      include 'strings.f'
      include 'splinecode.f'

c     eof

