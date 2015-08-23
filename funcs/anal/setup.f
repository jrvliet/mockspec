c     PACKAGE: Mockspec/anal
c     PROGRAM: sysanal
c     MODULE:  setup.f
c
c     DESCRIPTION 
c     routines that do the all important book keeping 
c
c
c     this file contains:
c     SUBROUTINE gettransitions
c     SUBROUTINE setup
c     SUBROUTINE sortions
c
c
c.........................................................................
c

      SUBROUTINE gettransitions(iknt,tranilist,losnum,specname,
     &                  ntran,trani,transpecfile,losindex,error)

c     populate the transition list from the exisiting files and store
c     relelvanet survey information about the transition

c     this routine is admittedly a little messy for it tries to take
c     care of several book keeping steps in the fly during a file read;
c     so I would recommend not messing with it (or if you do, save a
c     copy first!)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'sysanal.h'

      logical            error
      integer            j,iflag,anum,stage,ntran(mxions),iknt
      integer            access
      double precision   w,f
      character*80       tranilist,testfile,losnum,sname,ionname,
     &                   ionlabel,traniname,fileroot,old_iname,header,
     &                   specname(mxions)
      character*80       losindex,trani(mxions,mxtrans),
     &                   transpecfile(mxions,mxtrans)

      include            'sysanal.com'


c     basically, we cull through the transition list on the fly and
c     check for desired transition using the "iflag=1" condition; we
c     then construct the "QS0.*.gas_lines" file names corresponding to
c     the transition (i.e., "HI" for transition "Lya", "HI1025",
c     etc). these files contain the absorption line data to be modeled;
c     we then check if the file has been created (on disk); if they have
c     then we want to make the spectra for the desired transition of
c     this species so we increment the transition counter (ntran) and
c     store the name of the transition (trani) and the name of the file
c     the line data are stored in (linelist)


      error = .false.

c     IKNT is hard set to = 1, since we do only one species at a time on
c     the fly

      iknt  = 1

      DO 07 i=1,mxions
       ntran(i)  = 0
 07   CONTINUE


      CALL substr1(losnum,ionlabel)

      OPEN(unit=44,file=tranilist,status='old')  ! open master transition list
      READ(44,*,ERR=99) header                   ! stip header

      DO 21 j=1,mxions*mxtrans                   ! loop and read entries
       READ(44,*,END=22,ERR=99) iflag,sname,anum,stage,ionname,
     &                          traniname,w,f  

       IF ((ionname.eq.ionlabel).AND.(iflag.eq.1)) then  ! do this one
        CALL fappend2(losnum,traniname,'spec',testfile)  ! make .spec file name
        IF (access(testfile,'r').eq.0) then              ! file exist? then book keep
         specname(iknt)            = ionname
         ntran(iknt)               = ntran(iknt) + 1
         trani(iknt,ntran(iknt))   = traniname
         wave0(iknt,ntran(iknt))   = w
         fosc(iknt,ntran(iknt))    = f
         transpecfile(iknt,ntran(iknt)) = testfile
         DO 19 i=1,nions
          IF ((sname.eq.element(i)).AND.(stage.eq.ionstage(i))) then
           instrid(iknt) = instr(i)
           EWcut(iknt)   = EWlim(i)
           Slevid(iknt)  = Nsigion(i)
          END IF
 19      CONTINUE
        ENDIF
       END IF

 21   CONTINUE  ! read the next entry in the master transition list

 22   CLOSE(unit=44)


c     perform a sum check; if fail communiate and terminate

C     IF (ntran(iknt).lt.2) then
C      error = .true.
C      WRITE(1,602) losindex,specname(iknt)
C      WRITE(6,602) losindex,specname(iknt)
C     END IF

      RETURN

c     formatting error?

 99   WRITE(6,*) 'ERROR(gettransitions): atomic/transition data file ',
     &            tranilist(1:40)
      WRITE(6,*) 'is apparently not formatted correctly. Please'
      WRITE(6,*) 'read the Mockspec/Mockspec.help file for details.'
      CLOSE(unit=44)
      STOP

 602  FORMAT(1x,a8,2x,a9,2x,'-> no spectrum for this los. skipping.')

      END


c
c.........................................................................
c

      SUBROUTINE setup(ioni,ntran,trani,transpecfile,instrlist)

c     this routine first translates all the input into the required
c     architecture to run sysanal 9based upon original architecture for
c     more general cases), then it ties transitions together, and
c     finally, it gets the instrumental parameters

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i,j,k,ioni,ntran(mxions)
      double precision   R_isf
      character*80       spectrograph,instrlist,header
      character*80       trani(mxions,mxtrans),
     &                   transpecfile(mxions,mxtrans)

      include            'sysanal.com'



c     set up the arrays according to the architecture (left over from
c     the orginal sysanal)

      norders = ntran(ioni)

      DO 21 j=1,norders                 
       tie(j)       = 1                 ! all ions tied with tie=1
       N_sigma(j)   = Slevid(ioni)      ! minimum value
       order(j)     = trani(ioni,j) 
       specfiles(j) = transpecfile(ioni,j) 
       lambda0(j)   = wave0(ioni,j)
       fval(j)      = fosc(ioni,j)
 21   CONTINUE

c     architecture requires that the first tranisition in the arrays is
c     the strongest transition; so sort the transitions strength

      ntie = 1   ! hard wire all tranistions tied

      IF (norders.gt.1) then 

        CALL sortions(Slevid(ioni))

c     set up the transition ties; ntie is the number of unique species
c     (i.e., MgII, FeII, etc); tie(order) is the tie index of the tie
c     for that species, and njtie is the number of ion transitions for
c     the tie index; set the tie information

        DO 43 i=1,ntie
          njtie(i) = 0
         DO 41 j=1,norders
           IF (tie(j).eq.i) njtie(i) = njtie(i) + 1
 41      CONTINUE
 43     CONTINUE

        DO 47 i=1,ntie
          k = 0
          DO 45 j=1,norders
            IF (tie(j).eq.i) THEN
              k = k + 1
              jtie(i,k) = j
            END IF
 45       CONTINUE     
 47     CONTINUE

      ELSE

        jtie(1,1) = 1
        njtie(1)  = 1

      ENDIF

c     set up the instrument paramters for this ion; basically, we need
c     the resolution

      OPEN(unit=44,file=instrlist,status='old')

      READ(44,*) header
      DO 11 i=1,100
        READ(44,*,END=12,ERR=99) spectrograph,R_isf
        IF (spectrograph.eq.instrid(ioni)) then
         profile  = 1.0d0 / (2.35d0*R_isf)
         GOTO 12
        END IF
 11   CONTINUE

 12   CLOSE(unit=44)

c     return

      RETURN

c     formatting error?

 99   WRITE(6,*) 'ERROR(setup): instrument configuration file ',
     &            instrlist(1:40)
      WRITE(6,*) 'is apparently not formatted correctly. Please'
      WRITE(6,*) 'read the Mockspec/Mockspec.help file for details.'
      CLOSE(unit=44)
      STOP

      END

c
c.........................................................................
c

      SUBROUTINE sortions(Slevion)

c     the transition list can be a long table and there are no
c     constrains on the order in which transition were listed for a
c     given species; the architecture of sysanal requires that the first
c     ion in the list is the strongest transition, a good example would
c     be the P_3/2 - P_1/2 member of a resonant doublet as compared to
c     the P_1/2 - P_1/2 member.  the relative strengths are given by
c     f*lambda

c     here we find the index order in which the value f*lambda is
c     ascending and then apply this in reverse index order to sort the
c     accompanying relevant arrays of the ions in descending order

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i,j,k,nlist,indx(mxions)
      double precision   tmpNsig(mxions),tmplam(mxions),tmpf(mxions)
      double precision   flam(mxions),Slevion
      character*80       tmpord(mxions),tmpsfil(mxions)
      include            'sysanal.com'



c     local variable of the number of ions

      nlist = norders

c     store the arrays in temporary arrays to avoid overwriting; while
c     stroing, compute the quantity f*lambda for the sorting order

      DO 09 j=1,nlist
       tmpNsig(j) = N_sigma(j)   
       tmpord(j)  = order(j)     
       tmpsfil(j) = specfiles(j) 
       tmplam(j)  = lambda0(j)   
       tmpf(j)    = fval(j)      
       flam(j)    = fval(j)*lambda0(j)
 09   CONTINUE

c     get the index ordering of the flam array; returned in ascending
c     order

      CALL indexx(nlist,flam,indx)

c     sort in descending order using index k to index the indx array
c     in reverse order

      DO 11 j=1,nlist 
       k = nlist - (j-1)          ! reverse the indesing order
       i = indx(k)
       N_sigma(j)   = tmpNsig(i)  
       order(j)     = tmpord(i)   
       specfiles(j) = tmpsfil(i)  
       lambda0(j)   = tmplam(i)   
       fval(j)      = tmpf(i)     
 11   CONTINUE

c     now assign the survey significance limit to the strongest
c     transitions

      N_sigma(1) = Slevion
      
c     debug check

C      DO 13 j=1,nlist 
C       WRITE(6,*) N_sigma(j),order(j)(1:10),specfiles(j)(1:30),
C     &            lambda0(j),
C     &            fval(j)     
C 13   CONTINUE



      RETURN
      END
