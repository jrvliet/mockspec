c
c.........................................................................
c

      SUBROUTINE outfiles(kdo,jdo,GZoutfile,GZcoolfile)

c     this routine reads "rates.outfiles" and sets up which ions are to
c     be output; also sets up the cooling time data file for output

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      integer            i,ipos,knt
      integer            kdo(Imax),jdo(Imax)
      character*80       ionname,header,GZcoolfile,GZoutfile(Imax)
      include           'rates.com'

c     file "rates.outfiles" must be a local file

      OPEN(unit=17,file='rates.outfiles',status='old')
      DO i=1,noutfiles
       READ(17,*,END=99,ERR=99) GZoutfile(i),kdo(i),jdo(i)
      ENDDO
      CLOSE(unit=17)

c     construct the cooling time data file from the fist entry of the
c     GZoutile array; we append "tcdat.txt" following the 2nd "."

      ipos = 0
      knt  = 0
      DO i=1,80
       IF (GZoutfile(1)(i:i).eq.'.') then
        knt = knt + 1
        IF (knt.eq.2) then
         ipos = i
         GOTO 01
        ENDIF
       ENDIF      
      ENDDO

 01   GZcoolfile = GZoutfile(1)(1:ipos)//'tcdat.txt'
   
      RETURN

 99   WRITE(6,*) 'ERROR(outfiles): trouble reading file'
      WRITE(6,*) 'rates.outfiles.  did you correctly specify'
      WRITE(6,*) 'the number of output files in rates.inp,'
      WRITE(6,*) 'or is tere a formatting problem?'
      STOP

      END

c
c.........................................................................
c

      SUBROUTINE readionIDs

c     this routine reads the file 'tab-ionkj-ID.dat' and stores the
c     ionID array indexed by k,j; ex: ionID(k,j) for k=12 j=2 is
c     ionID(12,2)='MgII'

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      integer            i,j,k,kk,m,nel
      character*80       ionname,header
      character*250      fpath,infile
      include           'rates.com'


c     zero the vectorized index counter

      m = 0

c     ............
c     READ ION IDs
c     ............

c     create the table file name
      fpath = tabpath
      CALL sappend(fpath,'tab-ionkj-ID.dat',infile)

      OPEN(unit=1,file=infile,ERR=998,status='old')
      READ(1,*) header

      DO 11 i=1,Nmax
       READ(1,*,END=12) k,nel,j,ionname
       ionID(k,j) = ionname
       IF (j.eq.1.AND.specflag(k)) then
        m       = m + 1
        kidx(m) = k
       END IF
 11   CONTINUE

 12   CLOSE(unit=1)

      IF (ionname.ne.'ZnXXXI') then
       WRITE(6,*) 'ERROR(readionIDs): end of file = ionkj,dat'
       WRITE(6,*) 'not reached; should be Nmax-1 entries'
       WRITE(6,*) 'last entry read: ',i,k,j,ionID(k,j)(1:10)
       STOP
      END IF

c     define the number of species being modeled

      Nspecies = m


c     ................
c     READ SPECIES IDs
c     ................

c     create the table file name
      CALL sappend(fpath,'tab-speck-ID.dat',infile)

      OPEN(unit=1,file=infile,ERR=999,status='old')
      READ(1,*) header

      DO 31 k=1,Imax
       READ(1,*,END=32) kk,ionname
       IF (specflag(k)) specID(k) = ionname
 31   CONTINUE

 32   CLOSE(unit=1)


      RETURN


 998  WRITE(6,*) 'ERROR(readionIDs): Ion indexing file not found'
      WRITE(6,*) infile
      STOP

 999  WRITE(6,*) 'ERROR(readionIDs): Species indexing file not found'
      WRITE(6,*) infile
      STOP

      END


c
c.........................................................................
c

      SUBROUTINE readshellIDs

c     this routine reads the file 'tab-shell-ID.dat' and stores the
c     shell ID indexed by the shell number (see comments in main
c     routine)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      integer            s,i,n,l
      character*80       header
      character*250      fpath,infile
      include           'rates.com'


c     create the table file name
      fpath = tabpath
      CALL sappend(fpath,'tab-shell-ID.dat',infile)

      OPEN(unit=1,file=infile,ERR=999,status='old')
      READ(1,*) header
      DO 11 s=1,NSHmax
       READ(1,*) i,n,l,shellID(s)
 11   CONTINUE

      CLOSE(unit=1)

      RETURN

 999  WRITE(6,*) 'ERROR(readshellIDs): shell indexing file not found'
      WRITE(6,*) infile
      STOP

      END


c
c.........................................................................
c

      SUBROUTINE readeconfigs

c     this routine reads the file 'tab-isoseq-ID.dat' and stores the
c     ions isoelectronic sequence, Group on the periodic table, and the
c     full electron configuration; example

c     BeII: iso-sequence = 'lithium', Group = 'IA' and electron
c     configuration = '1s2.2s1'

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      integer            i,nel
      character*80       header
      character*250      fpath,infile
      include           'rates.com'

c     create the table file name
      fpath = tabpath
      CALL sappend(fpath,'tab-isoseq-ID.dat',infile)

      OPEN(unit=1,file=infile,ERR=999,status='old')
      READ(1,*) header
      DO 11 nel=1,Imax
       READ(1,*) i,sequence(i),group(i),config(i)
 11   CONTINUE

      CLOSE(unit=1)

      RETURN

 999  WRITE(6,*) 'ERROR(readeconfigs): isoelectronic file not found'
      WRITE(6,*) infile
      STOP

      END


c
c.........................................................................
c

      SUBROUTINE openfiles(GZoutfile,GZcoolfile)

c     this opens the output files 
c
c     note that the unit numbers are fixed, so that the WRITE and CLOSE
c     statements to these files must match the unit numbers used here
c     (the WRITE statements are hard wired in the subroutines in the
c     module "getcube.f"; the CLOSE statements are in the subroutine
c     "cfiles" below this subroutine)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      integer            j
      character*80       GZcoolfile,GZoutfile(Imax) 
      include           'rates.com'

c     open the ion files 

      DO 91 j=1,noutfiles
       OPEN(unit=j+10,file=GZoutfile(j),status='unknown')
 91   CONTINUE

c     open the cooling data file

      OPEN(unit=77,file=GZcoolfile,status='unknown')

      RETURN
      END


c
c.........................................................................
c

      SUBROUTINE closefiles

c     this closes the output files 
c
c     note that the unit numbers are fixed, so that the OPEN and WRITE
c     statements to these files must match the unit numbers used here
c     (this is hard wired in the subroutines in the module "getcube.f")

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      integer            j
      include            'rates.com'

c     close the ion files 

      DO 91 j=1,noutfiles
       CLOSE(unit=j+10)
 91   CONTINUE

c     close the cooling data file

      CLOSE(unit=77)

      RETURN
      END
