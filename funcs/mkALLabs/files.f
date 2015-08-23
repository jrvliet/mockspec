c     PACKAGE: Mockspec
c     PROGRAM: cullabs
c     MODULE:  files.f
c
c     this file contains:
c     SUBROUTINE ckfiles
c     SUBROUTINE mkfiles
c     SUBROUTINE openALLfiles
c     SUBROUTINE wrtdata
c     SUBROUTINE cullabsdata
c     SUBROUTINE cullregdata
c
c
c.............................................................................
c

      SUBROUTINE ckfiles(qsolist,tranilist)

c     check that the files entered on the command line can be accessed
c     in the present working directory; this routine does not confirm
c     the formating of the contenst of the file; only that they can be
c     accessed; for formatting contraints, see the file
c     "Mockspec/Mockspec.help"

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit          none
      logical           error
      integer           istat,access
      character*80      qsolist,paramlist,tranilist


      error = .false.

c     checking command line entry $1

      istat = access(qsolist,'r')
      IF (istat.ne.0) then
       WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &             qsolist(1:40)
       error = .true.
      END IF

c     checking command line entry $2

      istat = access(tranilist,'r')
      IF (istat.ne.0) then
       WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &             tranilist(1:40)
       error = .true.
      END IF

c     terminate if an error occured

c       error trap- terminate program

      IF (error) then
       WRITE(6,*) '***********************************************'
       WRITE(6,*) '* 1 or more of the command line files were not *'
       WRITE(6,*) '* present in the present working directory.    * '
       WRITE(6,*) '* For help, read the "Mockspec/Mockspec.help"  *'
       WRITE(6,*) '************************************************'
       STOP
      END IF

c     happy! return and GO!

      RETURN

      END



c
c..............................................................................
c
      SUBROUTINE mkfiles(klos)

c     populate the sysabsfile and regabsfile names
c     these are the files to be READ in

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'cullabs.h'

      integer            klos
      character*80       infile,outfile

      include            'cullabs.com'

      infile = losfile(klos)

      CALL fappend1(infile,'sysabs',outfile)
      sysabsfile(klos) = outfile
 
      CALL fappend1(infile,'regabs',outfile)
      regabsfile(klos) = outfile

      RETURN

      END


c
c......................................................................
c

      SUBROUTINE closeALLfiles

c     this closes the sysabs and regabs files

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include           'cullabs.h'

      integer           i,icol,funit
      character*80      infile,outfile1,outfile2

      include           'cullabs.com'

      CLOSE(unit=41)
      CLOSE(unit=42)

      RETURN

      END


c
c......................................................................
c

      SUBROUTINE openALLfiles

c     this routine makes the name and opens the output files (ALL files)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include           'cullabs.h'

      integer           i,icol,funit
      character*80      infile,outfile1,outfile2,aval

      include           'cullabs.com'


c     grab the expansion parameter

      infile = '../gal_props.dat'
      OPEN(unit=87,file=infile,status='old')
       DO i=1,6
        READ(87,*) outfile1  ! dummy reads
       ENDDO
       READ(87,*) aval       ! the a value (dummy variable)
      CLOSE(unit=87)

      infile = losfile(1)

c     create the ALL sysabs file name and open

      CALL fappend3(infile,'ALL',aval,'sysabs',outfile1)
      OPEN(unit=41,file=outfile1,status='unknown')

c     create the ALL regabs file name and open

      IF (nsubsys.gt.0) then
       CALL fappend3(infile,'ALL',aval,'regabs',outfile2)
       OPEN(unit=42,file=outfile2,status='unknown')
      ENDIF

      RETURN

      END


c......................................................................
c

      SUBROUTINE wrtdata

c     this routine writes the sysabs data to the ALL sysabs files

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include           'cullabs.h'

      integer           icell,ireg,icol,i,j,k,m,funit

      include           'cullabs.com'



c     loop over the unique ions; grad the relevant indices from the
c     translation table and write the sysabs data to the ALL file; when
c     done, close the file

c     the *.ALL.sysabs files

      WRITE(41,1000) 'los','D',headline

      DO 21 k=1,nlos
       IF (nreg(k).ne.0) then
        DO 23 ireg=1,nreg(k)
         WRITE(41,1200) losID(k),impact(k),sysdata(k,ireg)
 23     CONTINUE
       END IF
 21   CONTINUE


c     the *.ALL.regabs files

      IF (nsubsys.gt.0) then 

       WRITE(42,1200) losID(k),impact(k),sysdata(k,ireg)

       DO 31 k=1,nlos
        IF (nsubreg(k).ne.0) then
         DO 33 ireg=1,nsubreg(k)
          WRITE(42,1200) losID(k),impact(k),regdata(k,ireg)
 33      CONTINUE
        END IF
 31    CONTINUE

      END IF


      RETURN

 1000 FORMAT(1x,a3,2x,a6,1x,a190)
 1200 FORMAT(1x,a4,1x,f6.1,1x,a190)
     
      END



c
c......................................................................
c

      SUBROUTINE cullabsdata(klos)

c     this routine is called one los at a time, so we only work on one
c     los and then return to main; losnum is incremented in main

c     READ in the sysabs data

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include           'cullabs.h'

      integer           ireg,klos
      character*250     header,stringline,dummyline

      include           'cullabs.com'


      dummyline = '  0.0000000      0.00      0.00    0.000    0.000'//
     &            '    0.000    0.000    0.000      0.00      0.00'//
     &            '      0.00      0.00      0.00      0.00    0.000'//
     &            '    0.000    0.000    0.000    0.000    0.000'

c     set the maximum allowed number of columns to a local variable

      OPEN(unit=4,file=sysabsfile(klos),err=100,status='old')
      READ(4,'(a250)') header    ! header
      DO 09 ireg=1,mxreg
       READ(4,'(a250)',END=200) stringline     ! read data
       sysdata(klos,ireg) = stringline
 09   CONTINUE

      CLOSE(unit=4)               
      WRITE(6,*) 'ERROR(cullabsdata): file truncated; increase mxreg'
      STOP

 200  CLOSE(unit=4)
      nreg(klos) = ireg - 1
      IF (ireg.gt.0) headline = header  ! apparently not used
      RETURN

 100  CLOSE(unit=4)
      nreg(klos) = 1
      sysdata(klos,1) = dummyline
      RETURN

      END



c
c......................................................................
c

      SUBROUTINE cullregdata(klos)

c     this routine is called one los at a time, so we only work on one
c     los and then return to main; losnum is incremented in main

c     READ in the regabs data
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include           'cullabs.h'

      integer           ireg,klos
      character*250     stringline

      include           'cullabs.com'


c     set the maximum allowed number of columns to a local variable

      OPEN(unit=4,file=regabsfile(klos),err=100,status='old')
      READ(4,'(a250)') stringline    ! header
      DO 09 ireg=1,mxreg
       READ(4,'(a250)',END=200) stringline     ! read data
       regdata(klos,ireg) = stringline
 09   CONTINUE

      CLOSE(unit=4)                     
      WRITE(6,*) 'ERROR(cullregdata): file truncated; increase mxreg'
      STOP

 200  CLOSE(unit=4)
      nsubreg(klos) = ireg - 1

 100  RETURN

      END

