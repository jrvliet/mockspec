c     PACKAGE: Mockspec
c     PROGRAM: cullabs
c     MODULE:  main driver
c
c     INCLUDED MODULES:
c     access.f
c     files.f
c     subs.f
c     parze.f
c     index.f
c
c
c     DESCRIPTION 
c     this program takes as input the absorption line data from the
c     outputs of the program 'Mockspec/anal/sysanal' and culls together
c     the absorption line data and outputs master files that contain all
c     the information for all sightlines for each ionic species (i.e.,
c     HI, MgII, etc)
c
c     the code is very generalized in that no variables are specifically
c     assigned; it is strickly list directed and the variables are
c     symbolic; as such, if upstream modifications to the output files
c     from Mockspec/anal/sysanal are made then this program *should*
c     still work; howerve, it may be possible that parameters in the
c     cullabs.h file require modification if the number of columns are
c     increased (for example); sanity checks are performed where in key
c     locations to warn the user when this may be necessary
c
c
c     author: Chris Churchill
c     email:  cwc@nmsu.edu
c
c     last update: May 2011
c
c
c..............................................................................
c

      PROGRAM CULLABS

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include           'cullabs.h'

      integer           i,klos

      character*80      qsolist,tranilist
      character*80      losnum

      include           'cullabs.com'




c     grab the command line arguments

      CALL getarg(1,qsolist)     ! list of sightlines (QS0.XXX.dat files)

      IF (qsolist.eq.'') THEN
       WRITE(6,*) 'Usage: cullabs $1 [$2] '
       WRITE(6,*) '  $1  = filename of list of QS0.XXX.dat files'
       WRITE(6,*) ' [$2] = atomic/transition data file'
       WRITE(6,*) ' '
       WRITE(6,*) ' if $2 is null "Mockspec.transitions" is assumed.'
       WRITE(6,*) ' it is assumed this file has 1 header line.'
       WRITE(6,*) ' '
       WRITE(6,*) ' these files need to be in the present working'
       WRITE(6,*) ' directory: for required file formating, see'
       WRITE(6,*) ' the help file "Mockspec/Mockspec.help"'
       WRITE(6,*) ' '
       STOP
      END IF

      CALL getarg(2,tranilist)
      If (tranilist.eq.'') tranilist = 'Mockspec.transitions'

c     check that the qsolist and the tranilist are accessable in the
c     present working directory; if not communicate and bail

      CALL ckfiles(qsolist,tranilist)

c     clean up

      CALL clean


c
c     START
c

c     PART 1
c     loop over the qsolist; store the losnum file name for downstream
c     work; grab the impact parameter (global) of losnum in routine
c     getlosinfo; find, make and store the .sysabs and .regabs file
c     names (global) for this losnum in routine mkfiles, which also
c     returns the counter nions (global); store the header and system
c     absorption line data (global arrays) using routines cullabsdata
c     and cullregdata for the losnum; after culling through all losnum
c     close the qsolist file and set the nlos counter

c     initialize nions, the number of QS0.*:ion.sysabs and/or
c     QS0.*:ion.regabs files; nions is a global variable; intitialize
c     zero all the sysabs and regabs data matrices

      nsys    = 0
      nsubsys = 0

      CALL zeroall

      OPEN(unit=1,file=qsolist,status='old')
      OPEN(unit=45,file='Mockspec.runlog.cullabs',status='unknown')
      WRITE(6,600)
      WRITE(45,600)     

      DO 21 klos=1,maxlos
       nreg(klos)    = 0
       nsubreg(klos) = 0
       READ(1,*,END=22) losfile(klos)
       CALL getlosinfo(klos)
       CALL mkfiles(klos)
       CALL cullabsdata(klos)
       CALL cullregdata(klos)
 21   CONTINUE

c     error trap; if we get here then klos>maxlos; communicate and
c     terminate

      WRITE(6,*) 'ERROR(cullabs): The number of QSO sightlines in'
      WRITE(6,*) 'your qsolist exceeds the maximum number currently'
      WRITE(6,*) 'configured for program cullabs.  You can rectify'
      WRITE(6,*) 'the problem by editing cullabs.h and increasing'
      WRITE(6,*) 'the parameter MAXLOS; then re-make and re-run the'
      WRITE(6,*) 'cullabs program.'
      STOP

 22   CLOSE(unit=1)

      nlos = klos - 1

      DO i=1,nlos
       nsys    = nsys + nreg(i)
       nsubsys = nsubsys + nsubreg(i)      
       WRITE(6,602) losID(i),nreg(i),nsubreg(i)
       WRITE(45,602) losID(i),nreg(i),nsubreg(i)
      ENDDO


      WRITE(6,601)
      WRITE(45,601)     

c     PART 2
c     make the translation table so that we know which data get written
c     to which ALL:ion file file as indexed by the los index (klos) and
c     the ion index (j) in the below routines; open the ALL files and
c     write the data to the appropriate files

c     set the number of sightlines

      IF (nsys.gt.0) then
       CALL openALLfiles
       CALL wrtdata
       CALL closeALLfiles
      END IF

c     we are done; close the runlog file

      CLOSE(unit=45)

      STOP

 600  FORMAT(1x,'LOS',10x,'Nreg',2x,'Nsubreg',/,
     &       1x,'-----------',1x,'-----',2x,'------')
 601  FORMAT(1x,'-----------',1x,'-----',2x,'------')
 602  FORMAT(1x,a4,9x,i4,3x,i5)

      END


      include 'access.f'
      include 'files.f'
      include 'subs.f'
      include 'parze.f'
      include 'index.f'
