c     PACKAGE: Mockspec
c     PROGRAM: specsynth
c     MODULE:  files.f
c
c
c     DESCRIPTION
c     this module contains several subroutines that are involced in
c     accessing, reading, and storing data from files; several
c     routines are also dedicated to communication to the screen and
c     log files
c
c
c     this file contains:
c     SUBROUTINE ckfiles
c     SUBROUTINE getparams
c     SUBROUTINE getzabs
c     SUBROUTINE configspec
c     SUBROUTINE gettransitions
c     SUBROUTINE readlines
c     SUBROUTINE getHIseries
c     SUBROUTINE output
c     SUBROUTINE comm1
c     SUBROUTINE comm2
c     SUBROUTINE comm0lines
c     SUBROUTINE clean
c
c
c.............................................................................
c

        SUBROUTINE ckfiles(qsolist,paramlist,tranilist,instrlist)

c       check that the files entered on the command line can be accessed
c       in the present working directory; this routine does not confirm
c       the formating of the contenst of the file; only that they can be
c       accessed; for formatting contraints, see the file
c       "Mockspec/Mockspec.help"

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

	include           'specsynth.h'

        logical           error
        integer           istat,access
        character*80      qsolist,paramlist,tranilist,instrlist


        error = .false.

c       checking command line entry $1

        istat = access(qsolist,'r')

        IF (istat.ne.0) then
         WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &               qsolist(1:40)
         error = .true.
        END IF

c       checking command line entry $2

        istat = access(paramlist,'r')

        IF (istat.ne.0) then
         WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &               paramlist(1:40)
         error = .true.
        END IF

c       checking command line entry $3

        istat = access(tranilist,'r')

        IF (istat.ne.0) then
         WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &               tranilist(1:40)
         error = .true.
        END IF

c       checking command line entry $4

        istat = access(instrlist,'r')

        IF (istat.ne.0) then
         WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &               instrlist(1:40)
         error = .true.
        END IF

c       checking for the runlog of Mocksepc/mklos/los

        istat = access('Mockspec.runlog.los7','r')

        IF (istat.ne.0) then
         WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &              'Mockspec.runlog.los7 '
         WRITE(6,*) 'Have you run Mocksepc/mklos/los7'

         error = .true.
        END IF

c       error trap- terminate program

        IF (error) then
         WRITE(6,*) '***********************************************'
         WRITE(6,*) '* 1 or more of the command line files were not *'
         WRITE(6,*) '* present in the present working directory.    * '
         WRITE(6,*) '* For help, read the "Mockspec/Mockspec.help"  *'
         WRITE(6,*) '************************************************'
         STOP
        END IF

c       happy! return and GO!
        
        RETURN
        END



c.............................................................................
c

	SUBROUTINE getparams(paramlist)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

	include           'specsynth.h'

        integer           i
        double precision  dum
        character*80      paramlist,header



c       inititialize the number of ions we are working with

        nions = 0       

c       open the file (default="Mockspec.runpars") ; it is assumed to
c       have a single header line

	OPEN(unit=3,file=paramlist,status='old')
        READ(3,*) header
        DO 11, i=1,maxions
          READ(3,*,END=12,ERR=99) element(i),ionstage(i),dum,instr(i),
     &                            dum,snrat(i),vmax(i)
          nions = nions + 1
 11	CONTINUE

 12	CLOSE(unit=3)

c       required parameters, but not interactive anymore

        convolving = .false.
        slit     = 1.0
        resfac   = 3.0
        conwindo = 3

	RETURN


c       formatting error?

 99	WRITE(6,*) 'ERROR(getparams): survey configuration file ',
     &              paramlist(1:40)
        WRITE(6,*) 'is apparently not formatted correctly. Please'
        WRITE(6,*) 'read the Mockspec/Mockspec.help file for details.'
        CLOSE(unit=3)     
        STOP 

	END


c
c..............................................................................
c
      SUBROUTINE configspec(instrlist,j,wcen)

c     this routine grabs the instrumental parameters for this ion

c     the instrumental parameters are tabulated in the file
c     'Table.instruments' in the directory "/Mockspec/Codes/mkspec/"

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h'
      include           'const.dek'

      integer            i,j,k
      double precision   bigR,presel,readnoise,wcen
      character*80       instrlist
      character*80       header,spectrograph
      

c     we use "K" to denote the element ID, the index of the instrument
c     confiuration parameter

      k = elemid(j)

c     open, skip non-needed entries, and grab the galaxy redshift (zabs)

      OPEN(unit=44,file=instrlist,status='old')

      READ(44,*) header
      DO 11 i=1,100
        READ(44,*,END=12,ERR=99) spectrograph,bigR,presel,readnoise
        IF (spectrograph.eq.instr(k)) then
         vel_min  = -vmax(k) 
         vel_max  =  vmax(k)
         R_fac    = bigR
         rn       = readnoise
         snr      = snrat(k)
         GOTO 12
        END IF
 11   CONTINUE

 12   CLOSE(unit=44)     

c     compute the pixelization of the instrument

      wave_min = (1.0d0+vel_min/ckms)*(1.0d0+zabs)*lambda0(1)
      wave_max = (1.0d0+vel_max/ckms)*(1.0d0+zabs)*lambda0(1)
      wcen     = 0.50d0*(wave_max+wave_min) 
      dwave    = wcen/(presel*bigR)

      RETURN

c     formatting error?

 99   WRITE(6,*) 'ERROR(configspec): instrument configuration file ',
     &            instrlist(1:40)
      WRITE(6,*) 'is apparently not formatted correctly. Please'
      WRITE(6,*) 'read the Mockspec/Mockspec.help file for details.'
      CLOSE(unit=44)     
      STOP 

      END


c
c..............................................................................
c

      SUBROUTINE gettransitions(tranilist,losnum,ntran,trani,linelist)

c     populate the transition list and atomic data

c     aprt from file matching; there is no sanity check here for
c     transitions matching the ionic species; since SPECSYNTH is run
c     following the LOS program, and SPECYNTH used the output of LOS,
c     there should be no clash; in fact, the sanity check that at least
c     one transition will be set high for each ionic species in the
c     "transitions.dat" file is performed in the LOS program

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h'

      integer            i,j,iflag,anum,stage,ntran
      integer            access
      double precision   w0,fosc,damp
      character*80       tranilist,header,ionlabel
      character*80       testfile,losnum,sname,ionname,traniname
      character*80       trani(maxtrans),linelist(maxtrans)
      

c     basically, we cull through the transition list on the fly and
c     check for desired transition using the "iflag=1" condition; we
c     then construct the "*.lines" file names corresponding to the
c     transition (i.e., "OVI" for transitions). these files contain the
c     absorption line data to be modeled; we then check if the file has
c     been created (on disk); if it has then we want to make the spectra
c     for the desired transition of this species so we increment the
c     transition counter (ntran) and store the relevant data for each
c     transition

      OPEN(unit=44,file=tranilist,status='old')
      READ(44,*,ERR=99) header

      DO 11 j=1,maxtrans
       READ(44,*,END=12,ERR=99) iflag,sname,anum,stage,ionname,
     &                          traniname,w0,fosc,damp
       IF (iflag.eq.1) then
        CALL substr1(losnum,ionlabel)
        CALL fappend1(losnum,'lines',testfile)
        IF ((access(testfile,'r').eq.0).AND.   ! file exist?
     &      (ionlabel.eq.ionname)) then        ! ion match?
         ntran           = ntran + 1   ! increment transition counter
         trani(ntran)    = traniname   ! transition name
         linelist(ntran) = testfile    ! *.lines file name 
         lam0(ntran)     = w0          ! trani central wavelength 
         f0(ntran)       = fosc        ! trani oscillator strength
         gamma0(ntran)   = damp        ! trani damping constant
         elemid(ntran)   = 0
        ELSE
         GOTO 11
        END IF 
        DO 09 i=1,nions
         IF ((sname.eq.element(i)).AND.(stage.eq.ionstage(i))) then
          elemid(ntran) = i
         END IF
 09	CONTINUE
       END IF
 11   CONTINUE
 12   CLOSE(unit=44)     

      IF (ntran.eq.0) then
       WRITE(6,*) 'ERROR(gettransitions): could not identify'
       WRITE(6,*) 'transitions'
       STOP
      ENDIF


      RETURN

c     formatting error?

 99   WRITE(6,*) 'ERROR(gettransitions): atomic/transition data file ',
     &            tranilist(1:40)
      WRITE(6,*) 'is apparently not formatted correctly. Please'
      WRITE(6,*) 'read the Mockspec/Mockspec.help file for details.'
      CLOSE(unit=44)     
      STOP 


      END


c
c..............................................................................
c
      SUBROUTINE readlines(linelist,traniname)

c     read designated linelist files, which we have already verified
c     are available on disk in routine gettransition

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c


      include           'specsynth.h'

      integer           k,istat
      integer           access
      character*80      linelist,header,traniname


c     we will not bother to check the file format, since the use
c     should in no way have modified these output files from
c     program Mockspec/mklos/los

      lines = 1

      OPEN(unit=2,file=linelist,status='old')

      READ(2,*) zabs

      DO 11 while (lines.le.maxlines)

       k = lines
       read(2,*,end=12) zline(k),nline(k),bline(k)
       ion_name(k) = traniname
       IF (nline(k).ge.0.0d0) then
          nline(k) = 10.0d0**(nline(k)-13.0d0)
       ELSE
          nline(k) = -10.0d0**(abs(nline(k))-13.0d0)
       END IF
       IF (INDEX(ion_name(k),'HIs').ne.0) then
         CALL getHIseries
         lines = lines + 32
       ELSE
         lines = lines + 1
       END IF

 11   CONTINUE

 99   WRITE(6,*) 'WARNING(readlines): the Mockspec/mklos/los outfile ',
     &            linelist(1:50)
      WRITE(6,*) 'has more entries then the max presently configured'
      WRITE(6,*) 'for Mockspec/mksepec/specsynth to handle.  This can'
      WRITE(6,*) 'be rectified by editing the specsynth.h file, then '
      WRITE(6,*) 'increasing the value of parameter MAXLINES, and then'
      WRITE(6,*) 'typing "make" to recompile specynth.'

 12   CLOSE(unit=2)

c     due to the while loop logic, we correct the line count

      lines = lines - 1

      RETURN

      END


c
c.............................................................................
c

      SUBROUTINE getHIseries

c     though there is no file I/I here, we do access a table (and
c     assumption are implied here) and this is called from a I/O
c     routine

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h'
      integer           i,j,k

      character*80      HIseries(32)


      data HIseries / 'HI1026', 'HI972' , 'HI950' , 'HI938' , 'HI931' , 
     @                'HI926' , 'HI923' , 'HI921' , 'HI919' , 'HI918' , 
     @                'HI917' , 'HI916a', 'HI916b', 'HI915a', 'HI915b',
     @                'HI915c', 'HI914a', 'HI914b', 'HI914c', 'HI914d',    
     @                'HI913a', 'HI913b', 'HI913c', 'HI913d', 'HI913e',    
     @                'HI913f', 'HI913g', 'HI913g', 'HI913h', 'HI913i',    
     @                'LLS'   , 'Lya'   
     @              /

      j = lines
      k = lines

      DO 11 i=1,32
        ion_name(k) = HIseries(i)
        zline(k)    = zline(j)
        nline(k)    = nline(j)
        bline(k)    = bline(j)
        k           = k + 1
 11   CONTINUE

      RETURN
      END


c
c..............................................................................
c

      SUBROUTINE output(losnum,trani)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h'
      include           'const.dek'
      integer           pixi
      double precision  wave,vel,flux,noise,w0z
      character*80      losnum,trani,delim,outfile


c     create the output file name

      CALL sappend(trani,'.spec',delim)
      CALL fappend1(losnum,'',outfile)
      CALL sappend(outfile,delim,outfile)

c     open the output file and write the spectrum

      OPEN(unit=3,file=outfile,status='unknown')

       DO 05 pixi=1,ndata
        wave  = lambda(pixi)
        w0z   = (1.0d0+zabs)*lambda0(1)
        vel   = ckms*(wave-w0z)/w0z
        flux  = wrkflx(pixi)
        noise = sigma(pixi)
        WRITE(3,100) wave,vel,flux,noise,noise,1.0d0
 05    CONTINUE

       CLOSE(unit=3)

      RETURN

 100  FORMAT(1x,f10.4,2x,f9.2,2x,1p4e13.5)

      END


c
c..............................................................................
c

      SUBROUTINE comm1(qsolist,paramlist,tranilist,instrlist)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      character*80       qsolist,paramlist
      character*80       tranilist,instrlist  


c     write to the runlog file

      WRITE(1,*) ' ' 
      WRITE(1,*) '*******************************************'
      WRITE(1,*) '- Mockspec/specsynth input files'
      WRITE(1,'(a,a40)') ' - LOS data list          : ',qsolist
      WRITE(1,'(a,a40)') ' - absline survey config  : ',paramlist
      WRITE(1,'(a,a40)') ' - atomic/transition data : ',tranilist
      WRITE(1,'(a,a40)') ' - instrument config data : ',instrlist
      WRITE(1,*) ' ' 
      WRITE(1,*) '*******************************************'

c     write to screen

      WRITE(6,*) ' ' 
      WRITE(6,*) '*******************************************'
      WRITE(6,*) '- Mockspec/specsynth input files'
      WRITE(6,'(a,a40)') ' - LOS data list          : ',qsolist
      WRITE(6,'(a,a40)') ' - absline survey config  : ',paramlist
      WRITE(6,'(a,a40)') ' - atomic/transition data : ',tranilist
      WRITE(6,'(a,a40)') ' - instrument config data : ',instrlist
      WRITE(6,*) ' ' 
      WRITE(6,*) '*******************************************'

      RETURN
      END


c
c..............................................................................
c

      SUBROUTINE comm2(losnum,j,trani,linelist)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h' 

      integer            j,k
      character*80       losnum,trani,linelist


      k = elemid(j)


c     if this is the first transition for this losnum, then write the
c     header


      IF (j.eq.1) then
       WRITE(1,*) ' '
       WRITE(1,600) 'LOS','tran','linelist','instr','R',
     &              'ang/pix','kms/pix','pix/res','kms/res','RN','S/N',
     &              'lam-','lam+','vel-','vel+','Npix','Ncells'
       WRITE(6,*) ' '
       WRITE(6,600) 'LOS','tran','linelist','instr','R',
     &              'ang/pix','kms/pix','pix/res','kms/res','RN','S/N',
     &              'lam-','lam+','vel-','vel+','Npix','Ncells'
      END IF


c     write to the runlog file

      WRITE(1,601) losnum,trani,linelist,instr(k),R_fac,
     &             dwave,dv,pixpres,profile,rn,snr,
     &             wave_min,wave_max,vel_min,vel_max,
     &             ndata,lines

c     write to screen

      WRITE(6,601) losnum,trani,linelist,instr(k),R_fac,
     &             dwave,dv,pixpres,profile,rn,snr,
     &             wave_min,wave_max,vel_min,vel_max,
     &             ndata,lines

      RETURN

 600  FORMAT(2x,a3,19x,a4,6x,a8,20x,a5,7x,a1,
     &       7x,a7,2x,a7,2x,a7,2x,a7,5x,a2,5x,a3,
     &       5x,a4,8x,a4,9x,a4,6x,a4,6x,a4,4x,a6)
 601  FORMAT(1x,a21,2x,a8,2x,a26,2x,a10,2x,f6.0,
     &       2x,f7.4,2x,f7.4,2x,f7.2,2x,f7.2,2x,f6.1,2x,f6.1,
     &       2x,f10.4,2x,f10.4,2x,f8.1,2x,f8.1,2x,i7,2x,i7)

      END

c
c..............................................................................
c

      SUBROUTINE comm0lines(losnum,j,trani,linelist)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h' 

      integer            j,k
      character*80       losnum,trani,linelist,note1,note2


      k = elemid(j)

      note1 = 'No cells found in linelist; see "Mockspec.los7.runlog"'
      note2 = 'and "Mockspec.los.errlog" files for details on bad cells'

c     if this is the first transition for this losnum, then write the
c     header


      IF (j.eq.1) then
       WRITE(1,*) ' '
       WRITE(1,600) 'LOS','tran','linelist','WARNING'
       WRITE(6,*) ' '
       WRITE(6,600) 'LOS','tran','linelist','WARNING'
      END IF


c     write to the runlog file

      WRITE(1,601) losnum,trani,linelist,note1,note2

c     write to screen

      WRITE(6,601) losnum,trani,linelist,note1,note2

      RETURN

 600  FORMAT(2x,a3,19x,a4,6x,a8,20x,a7)
 601  FORMAT(1x,a21,2x,a8,2x,a26,2x,a53,1x,a60)

      END

c
c.........................................................................
c

      SUBROUTINE clean

c     this routine deletes all files that will be written to; the
c     purpose is to avoid runtime errors that occur on some system
c     complaining that a file opened for writing cannot be written to
c     (sigh)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

c     remove the old runlog file if it exists

      CALL system('/bin/rm -f Mockspec.runlog.spec')

c     remove any old .spec files

      CALL system('/bin/rm -f QS0*.spec')

      RETURN
      END
