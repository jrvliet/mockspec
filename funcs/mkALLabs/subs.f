c     PACKAGE: Mockspec
c     PROGRAM: cullabs
c     MODULE:  subs.f
c
c     this file contains:
c     SUBROUTINE getlosinfo
c     SUBROUTINE fappend1
c     SUBROUTINE fappend2
c     SUBROUTINE getlosID
c     SUBROUTINE zeroall
c     SUBROUTINE clean
c
c
c......................................................................
c

      SUBROUTINE getlosinfo(klos)

c     get the expansion factor and impact parameter for the current los;
c     this routine is called one los at a time, so we only work on one
c     los and then return to main; klos and losnum are incremented in
c     main

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include           'cullabs.h'
      integer           n,klos
      double precision  field(maxcols),b1,b2
      character*250     stringline
      character*80      losnum,tmpstr
      include           'cullabs.com'


c     construct the file name and open

      n      = maxcols
      losnum = losfile(klos)

c     the second header line of all los files contains the sky locations
c     of the los; the x and y impact parameter positions; grab them

      OPEN(unit=2,file=losnum,status='old')
      READ(2,'(a250)') stringline  ! read first line and discard
      READ(2,'(a250)') stringline  ! read second line
      CALL parzeline(n,field,stringline)
      b1 = field(11)  ! impact on sky kpc (horizontal) is 6th field 
      b2 = field(12)  ! impact on sky kpc (vertical) is 8th field
      CLOSE(unit=2)

C     THIS WILL REQUIRE MODIFICATION 
C      impact(klos) = sqrt(b1*b1+b2*b2)
      impact(klos) = b1

c     now get the LOS ID number

      tmpstr = ' '
      CALL getlosID(losnum,tmpstr)
      losID(klos) = tmpstr

      RETURN

      END


c.........................................................................
c

      SUBROUTINE fappend1(infile,delim,outfile)

c     this routine takes infile, truncates off its delimeter and
c     replaces it with the two strings midstr and delim 

c     the delimiter of "infile" is assumed to be ".dat"; this is
c     truncated off the "infile" string and then two strings are
c     appended in the form ":midstrdelim"; for flexibility the "." is
c     not automatically placed between "midstr" and "delim"; the use
c     must include the "."  in the strong "delim"

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,k,lend
      character*(*)       infile,delim,outfile

c     lend returns the allocated string length, not the length of the
c     stored string

      lend = len(infile)
      k    = 0

      DO 09 i=1,lend
        k = i
        IF (infile(i:i+3).eq.'.dat') GOTO 10
 09   CONTINUE

 10   outfile = infile(1:k-1)//'.'//delim

      RETURN
      END

c.........................................................................
c

      SUBROUTINE fappend3(infile,midstr,aval,delim,outfile)

c     this routine takes infile, truncates off its delimeter and
c     replaces it with the two strings midstr and delim 

c     the delimiter of "infile" is assumed to be ".dat"; this is
c     truncated off the "infile" string and then two strings are
c     appended in the form ":midstrdelim"; for flexibility the "." is
c     not automatically placed between "midstr" and "delim"; the use
c     must include the "."  in the strong "delim"

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,jmid,java,jinf,lend,pknt
      character*(*)       infile,midstr,aval,delim,outfile


c     obtain the length of midstr string

      lend = len(midstr)
      jmid = 0
      DO i=1,lend
        jmid = i
        IF (midstr(i:i).eq.' ') GOTO 07
      ENDDO

c     obtain the length of the aval string

 07   lend = len(aval)
      java = 0
      DO i=1,lend
        java = i
        IF (aval(i:i).eq.' ') GOTO 08
      ENDDO

c     find the 2nd decimal place in the infile strong

 08   lend = len(infile)
      jinf = 0
      pknt = 0
      DO i=1,lend
        IF (infile(i:i).eq.'.') then
         pknt = pknt + 1
         IF (pknt.eq.2) then
          jinf = i
          GOTO 10
         ENDIF
        ENDIF
      ENDDO

c     make the outfile string

C 10   outfile = infile(1:jinf)//midstr(1:jmid)//'.a'//
C     &   aval(1:java-1)//'.'//delim

 10   outfile = infile(1:jinf)//'a'//aval(1:java-1)//'.'//
     &          midstr(1:jmid)//'.'//delim

      RETURN
      END


c.........................................................................
c

      SUBROUTINE getlosID(losnum,losID)

c     this routine takes the input string losnum, which is assumed to be
c     a *.dat file and culls out the XXX, which is the lOS ID number; it
c     assumed that it resides between "." 2 and "." 3

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,k1,k2,k3,lend
      character*(*)       losnum,losID

c     lend returns the allocated string length, not the length of the
c     stored string


      lend = len(losnum)
      k1   = 0
      k2   = 0

c     find the first period

      DO 09 i=1,lend
        k1 = i
        IF (losnum(i:i).eq.'.') GOTO 10
 09   CONTINUE

c     find the second period

 10   DO 11 i=k1+1,lend
        k2 = i
        IF (losnum(i:i).eq.'.') GOTO 12
 11   CONTINUE

c     the losID is the string between the periods

 12   k3 = 0
      DO 13 i=k2+1,lend
        k3 = i
        IF (losnum(i:i).eq.'.') GOTO 14
 13   CONTINUE


 14   losID = losnum(k2+4:k3-1)

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE zeroall

c     this routine nulls all the matrices so that elements that are
c     not reference later are explicitely nulled

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include             'cullabs.h'

      integer             klos,j,icol,ireg

      include             'cullabs.com'

   

      DO 04 icol=1,maxcols
       DO 03 klos=1,maxlos
        DO 02 j=1,mxreg
         sysdata(klos,j) = ' '
         DO 01 ireg=1,mxreg
          regdata(klos,j) = ' '
 01      CONTINUE
 02     CONTINUE
 03    CONTINUE
 04   CONTINUE


      RETURN

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

      CALL system('/bin/rm -f Mockspec.runlog.cullabs')

c     remove any old .sysabs files

      CALL system('/bin/rm -f ALL*.sysabs')

c     remove any old .regabs files

      CALL system('/bin/rm -f ALL*.regabs')

      RETURN
      END


c
