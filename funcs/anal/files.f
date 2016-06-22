c     PACKAGE: Mockspec/anal
c     PROGRAM: sysanal
c     MODULE:  files.f
c
c     DESCRIPTION
c     routines that manage the I/O to files and screen
c
c
c     this file contains:
c     SUBROUTINE ckfiles
c     SUBROUTINE getparams
c     SUBROUTINE getdata
c     SUBROUTINE comm1
c     SUBROUTINE comm0
c     SUBROUTINE comm2
c     SUBROUTINE writefiles
c
c
c.............................................................................
c

      SUBROUTINE ckfiles(qsolist,paramlist,tranilist,instrlist)

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
      character*80      qsolist,paramlist,tranilist,instrlist


      error = .false.

c       checking command line entry $1

      istat = access(qsolist,'r')
      IF (istat.ne.0) then
      WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &               qsolist(1:40)
       error = .true.
      END IF

c     checking command line entry $2

      istat = access(paramlist,'r')

      IF (istat.ne.0) then
       WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &               paramlist(1:40)
       error = .true.
      END IF

c     checking command line entry $3

      istat = access(tranilist,'r')
      IF (istat.ne.0) then
       WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &               tranilist(1:40)
       error = .true.
      END IF

c     checking command line entry $4

      istat = access(instrlist,'r')
      IF (istat.ne.0) then
       WRITE(6,*) 'ERROR(ckfiles): no access to file ',
     &               instrlist(1:40)
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

c     happy! return and GO!

      RETURN

      END


c.............................................................................
c

      SUBROUTINE getparams(paramlist)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'sysanal.h'

      integer           i,iknt
      double precision  dum
      character*80      paramlist,header

      include           'sysanal.com'



c     inititialize the number of ions we are working with

      nions = 0

c     open the file (default="Mockspec.runpars") ; it is assumed to have a
c     single header line

      OPEN(unit=3,file=paramlist,status='old')
      READ(3,*) header
      DO 11, i=1,mxions
        READ(3,*,END=12,ERR=99) element(i),ionstage(i),dum,instr(i),
     &                          EWlim(i),dum,dum,Nsigion(i)
        nions = nions + 1
 11   CONTINUE

 12   CLOSE(unit=3)

c     sum check that the paramlist matches the number of ions generated
c     by the LOS generating program, where nions is determined in
c     routine gettransitions (celled before this one


      RETURN


c     formatting error?

 99   WRITE(6,*) 'ERROR(getparams): survey configuration file ',
     &            paramlist(1:40)
      WRITE(6,*) 'is apparently not formatted correctly. Please'
      WRITE(6,*) 'read the Mockspec/Mockspec.help file for details.'
      CLOSE(unit=3)
      STOP

      END

c
c.........................................................................
c

      SUBROUTINE getdata(specname,error)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit   none
  
      include            'sysanal.h'

      logical            error
      integer            i,j,jj,k,kk,nline,ndum
      integer            access,istat
      parameter          (ndum = 20)
      double precision   dum(ndum)
      character*80       read_file,specname
      character*200      stringline

      include            'sysanal.com'


c     set the error low - it is set high if the file is not found, or it
c     is set high if there are no pixels in the file

      error = .false.

c     loop through all the ions for this species and store their
c     spectra; j is the transition index and i is the pixel index

      j = 1

      DO 11 WHILE (j.le.norders)

       ndata(j)  = 0
       read_file = specfiles(j)

       OPEN(unit=3,file=read_file,err=18,status='old')

       DO 15 i=1,nmx
        READ(3,'(a200)',END=17) stringline
        CALL parzeline(nline,ndum,dum,stringline)
        wave(i,j)  = dum(1)
        vel(i,j)   = dum(2)
        flux(i,j)  = dum(3)
        sigma(i,j) = dum(4)  
        smsig(i,j) = dum(5)  
        cont(i,j)  = dum(6)
 15    CONTINUE

c     if you get here, then the data was truncated warn and continue

       IF (iprint) then
        WRITE(1,600) specname,read_file,vel(nmx,j)
        WRITE(6,600) specname,read_file,vel(nmx,j)
       ENDIF

 17    CLOSE(unit=3)

       ndata(j) = i - 1

       IF (ndata(j).eq.0) then
         error = .true.
         RETURN
       END IF

       GOTO 19  ! increment order counter and continue

 18    IF (iprint) then
        WRITE(1,601) specname,read_file
        WRITE(6,601) specname,read_file
       ENDIF
       error = .true.

 19    j = j + 1

 11   CONTINUE

c     good return

      RETURN

c     formats

 200  FORMAT(1x,6f11.4)

 600  FORMAT(1x,a9,2x,'WARNING(getdata): NMX exceeded for ',a24,
     &       1x,'- spectrum truncated at vel > ',f10.2,' km/s')
      
 601  FORMAT(1x,a9,2x,'WARNING(getdata): no spectrum for ',a24,
     &       1x,'- skipping this ion')
      
      END


c
c..............................................................................
c

      SUBROUTINE comm1(qsolist,paramlist,tranilist,instrlist)

c     communicate the command line input files being used 

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      character*80       qsolist,paramlist
      character*80       tranilist,instrlist


c     write to the runlog file

      WRITE(1,*) ' '
      WRITE(1,*) '*******************************************'
      WRITE(1,*) '- Mockspec/sysanal input files'
      WRITE(1,'(a,a40)') ' - LOS data list          : ',qsolist
      WRITE(1,'(a,a40)') ' - absline survey config  : ',paramlist
      WRITE(1,'(a,a40)') ' - atomic/transition data : ',tranilist
      WRITE(1,'(a,a40)') ' - instrument config data : ',instrlist
      WRITE(1,*) ' '
      WRITE(1,*) '*******************************************'

c     write to screen

      WRITE(6,*) ' '
      WRITE(6,*) '*******************************************'
      WRITE(6,*) '- Mockspec/sysanal input files'
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

      SUBROUTINE comm0(j,losindex,spcnam)

c     communicate non detection for this ion; j=iorder

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'sysanal.h'

      integer            j
      double precision   Wcut,Siglev
      character*80       spcnam,spectrog,losindex,note

      include            'sysanal.com'

c     grab this ions quantities for writing

      spectrog = instrid(j)
      Wcut     = EWcut(j)
      Siglev   = Slevid(j)

c     set the note

      note = 'non detection'

c     write to the runlog file

      WRITE(1,601) losindex,spcnam,spectrog,Wcut,nlines,Siglev,note,
     &             (order(i),i=1,norders)


c     write to screen

      WRITE(6,601) losindex,spcnam,spectrog,Wcut,nlines,Siglev,note,
     &             (order(i),i=1,norders)


      RETURN

 601  FORMAT(1x,a8,2x,a9,2x,a10,2x,f6.4,2x,i3,6x,' ...  ',
     &       2x,f5.2,2x,a18,2x,10a10)

      END

c
c..............................................................................
c

      SUBROUTINE comm2(j,losindex,spcnam,detect)

c     communicate detection for this ion; j=iorder; but account for
c     whether the EW is above the survey threshold, EWcut (which is
c     stored in logical "detect"

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'sysanal.h'

      logical            detect
      integer            j
      double precision   Wcut,Siglev
      character*80       spcnam,spectrog,losindex,note

      include           'sysanal.com'



c     grab this ions quantities for writing

      restew   = ewtot(1)/(1.0+zbar)
      spectrog = instrid(j)
      Wcut     = EWcut(j)
      Siglev   = Slevid(j)


c     determine message

      IF (detect) then
       note = 'above EWcut'
C       note = 'detection included'
      ELSE
       note = 'below EWcut'
C       note = 'detection excluded'
      END IF

c     write to the runlog file

      WRITE(1,602) losindex,spcnam,spectrog,Wcut,nlines,restew,
     &             Siglev,note,(order(i),i=1,norders)

c     write to screen

      WRITE(6,602) losindex,spcnam,spectrog,Wcut,nlines,restew,
     &             Siglev,note,(order(i),i=1,norders)

      RETURN

 602  FORMAT(1x,a8,2x,a9,2x,a10,2x,f6.4,2x,i3,5x,f7.4,2x,
     &       f5.2,2x,a18,2x,10a10)

      END


c.........................................................................
c

      SUBROUTINE writefiles(losnum)

c
c     write all the files for all ion transitions
c
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit   none
  
      include            'const.dek'
      include            'sysanal.h'

      integer            i,j,k,j2
      double precision   rw,drw,v1,v2,asym,dasym
      character*80       write_file,losnum,fileroot

      include            'sysanal.com'


c     set J=strong member, J2=weaker member

      j  = 1
      j2 = 2

c     create the root name of the write_files

      IF (nlines.gt.1) THEN


c     create and open the file for the region data

      CALL fappend1(losnum,'regabs',write_file)
      OPEN(unit=13,file=write_file,status='unknown')
      WRITE(13,1300)

c     if more than one region...
c     write the velocity moments region by region 

      DO 25 i=1,nlines
        rw  = ew(i,j)/(1.0d0+zbar)
        IF (ewsig(i,j).ne.-1.0) then            ! limit ?
          drw = ewsig(i,j)/(1.0d0+zbar)         ! uncertainty in RW
        ELSE                                    ! or 
          drw = -1.0d0                          ! flag as limit 
        END IF      
        asym  = vasym(i,j)/(ckms*rw/lambda0(j))
        dasym = sigvasym(i,j)/(ckms*rw/lambda0(j))
        v1    = vel(f_beg(i,j,1),j)             ! blue vel extreme
        v2    = vel(f_end(i,j,1),j)             ! red vel extreme
        WRITE(13,1301) i,zbar,v1,v2,
     &        rw,drw,dr(i,j2),drsig(i,j2),
     &        siglevel(i,j),vbar(i,j),sigvbar(i,j),
     &        vwidth(i,j),sigvwidth(i,j),asym,dasym,
     &        0.0,1.25
 25   CONTINUE
      CLOSE(unit=13)

      END IF

c     create and open the file for the overall system results

      CALL fappend1(losnum,'sysabs',write_file)
      OPEN(unit=14,file=write_file,status='unknown')
      WRITE(14,1400)  

c     write the follwoing quantities to the file

c     1. the redshift, and equivalent widths
c
c     2. the total velocity moments; i.e., the mean velocity, velocity
c     width, and the asymmetry; we convert the asymmetry to asym/vew,
c     where vew = c*ew/lambda:
c
c     3. total measured optical depths and AOD column densities for each
c     ion transition; the errors are not symmetric (because of the
c     natural log), so there is a down error (dd) and an up error (du)
c     for each of the totals: "sysanal.aod"

       rw  = ewtot(j)/(1.0d0+zbar)                 ! rest frame EW
       IF (ewsigtot(j).ne.-1.0) then               ! limit ?
          drw = ewsigtot(j)/(1.0d0+zbar)           ! uncertainty in RW
       ELSE                                        ! or 
          drw = -1.0d0                             ! flag as limit 
       END IF      
       asym  = vtotasym(j)/(ckms*rw/lambda0(j))    ! vel asymmetry
       dasym = sigvtotasym(j)/(ckms*rw/lambda0(j)) ! uncertainty in ASYM
       v1    = vel(f_beg(1,j,1),j)                 ! blue vel extreme
       v2    = vel(f_end(nlines,j,1),j)            ! red vel extreme


       ddtautot(j) = log10((tautot(j)-ddtautot(j))/tautot(j))
       dutautot(j) = log10((tautot(j)+dutautot(j))/tautot(j))
       ddtautot(j) = abs(ddtautot(j))
       tautot(j)   = log10(tautot(j))

c      If SL is too high, prints ***** in column -> Cap the value
       IF (sltot(j).gt.1e3) then
          sltot(j) = 1e3
       END IF

       WRITE(14,1401) zbar,v1,v2,
     &       rw,drw,drtot(j2),drsigtot(j2),
     &       sltot(j),vtotbar(j),sigvtotbar(j),
     &       vtotwidth(j),sigvtotwidth(j),asym,dasym,
     &       tautot(j),ddtautot(j),dutautot(j),
     &       coltot(j),ddcoltot(j),ducoltot(j),
     &       0.0,1.25

      CLOSE(unit=14)

c     we are done looping through the ions

c     debugging

C      OPEN(unit=11,file='ew_regions.dat',status='unknown')
C      DO 27 i=1,nlines
C        w1 = wave(f_beg(i,1,1),1)
C        w2 = wave(f_end(i,1,1),1)
C        w3 = wbar(i,1)
C        WRITE(11,1110) i,1.25,w1,w2,w3,lambda0(1)
C 27   CONTINUE
C      CLOSE(unit=11)

c     this is a clean return

      RETURN

c     formats

 1300 FORMAT(1x,t4,'reg',t10,'zabs',t24,'v-',t34,'v+',
     &       t42,'EW_r',t50,'dEW_r',
     &       t60,'DR',t69,'dDR',t78,'SL',t88,'Vbar',
     &       t97,'dVbar',t107,'Vsprd',t115,'dVsprd',
     &       t127,'Vasym',t136,'dVasym',
     &       t145,'dum',t149,'ytick')
 1301 FORMAT(1x,i5,f10.7,2f10.2,5f9.3,6f10.2,f6.1,f6.2)
 1400 FORMAT(1x,t5,'zabs',t19,'v-',t29,'v+',t37,'EW_r',
     &       t45,'dEW_r',t55,'DR',t64,'dDR',t73,'SL',
     &       t83,'Vbar',t92,'dVbar',t102,'Vsprd',t111,'dVsprd',
     &       t122,'Vasym',t131,'dVasym',t142,'lgt',t150,'dtau-',
     &       t159,'dtau+',t169,'logN',t176,'dNcol-',t185,'dNcol+',
     &       t194,'dum',t198,'ytick')
 1401 FORMAT(1x,f10.7,2f10.2,5f9.3,6f10.2,6f9.3,f6.1,f6.2)

      END

