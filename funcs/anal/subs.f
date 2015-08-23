c     PACKAGE: Mockspec/anal
c     PROGRAM: sysanal
c     MODULE:  subs.f
c
c     DESCRIPTION 
c     routines that manage string and perform other book keeping
c     function
c
c
c     this file contains:
c     SUBROUTINE zeroall
c     SUBROUTINE sortregions
c     SUBROUTINE parzeline
c     SUBROUTINE parzstring
c     SUBROUTINE fappend
c     SUBROUTINE substr1
c     SUBROUTINE substr2
c     SUBROUTINE fappend1
c     SUBROUTINE fappend2
c     SUBROUTINE sappend
c     SUBROUTINE clean
c
c
c.........................................................................
c

      SUBROUTINE zeroall

c
c     zero all working arrays
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            j,i

      include            'sysanal.com'


c     intialize all working arrays


      nlines = 0

c     loop over orders

      DO 61 i=1,mxord

       ndata(i)        = 0
       ewtot(i)        = 0.0d0
       ewsigtot(i)     = 0.0d0
       drtot(i)        = 0.0d0
       drsigtot(i)     = 0.0d0
       sltot(i)        = 0.0d0
       vtotbar(i)      = 0.0d0
       vtotwidth(i)    = 0.0d0
       vtotasym(i)     = 0.0d0
       sigvtotbar(i)   = 0.0d0
       sigvtotwidth(i) = 0.0d0
       sigvtotasym(i)  = 0.0d0
       tautot(i)       = 0.0d0
       dutautot(i)     = 0.0d0
       ddtautot(i)     = 0.0d0
       coltot(i)       = 0.0d0
       ducoltot(i)     = 0.0d0
       ddcoltot(i)     = 0.0d0

c     feature/order arrays

       DO 71 j=1,mxlin

        sf_flag(i,j)   = .false.
        nfind(i,j)     = 0
        f_beg(j,i,j)   = 0
        f_end(j,i,j)   = 0
        ew(j,i)        = 0.0d0
        ewsig(j,i)     = 0.0d0
        wbar(j,i)      = 0.0d0
        sigwbar(j,i)   = 0.0d0
        width(j,i)     = 0.0d0
        sigwidth(j,i)  = 0.0d0
        dr(j,i)        = 0.0d0
        drsig(j,i)     = 0.0d0
        siglevel(j,i)  = 0.0d0
        vbar(j,i)      = 0.0d0
        vwidth(j,i)    = 0.0d0
        vasym(j,i)     = 0.0d0
        sigvbar(j,i)   = 0.0d0
        sigvwidth(j,i) = 0.0d0
        sigvasym(j,i)  = 0.0d0

 71    CONTINUE

c     pixel/order arrays

       DO 51 j=1,nmx

        mask(j,i)     = 1
        collim(j,i)   = 0
        wave(j,i)     = 0.0d0
        vel(j,i)      = 0.0d0
        flux(j,i)     = 0.0d0
        sigma(j,i)    = 0.0d0
        smsig(j,i)    = 0.0d0
        cont(j,i)     = 0.0d0
        ewpix(j,i)    = 0.0d0
        ewsigpix(j,i) = 0.0d0
        tau(j,i)      = 0.0d0
        dutau(j,i)    = 0.0d0
        ddtau(j,i)    = 0.0d0
        col(j,i)      = 0.0d0
        ducol(j,i)    = 0.0d0
        ddcol(j,i)    = 0.0d0


 51    CONTINUE


 61   CONTINUE


c     return

      RETURN
      END


c.........................................................................
c

      SUBROUTINE sortregions

c
c     here we sort the regions just in case they are not in velocity
c     order; keep track of all book keeping (whew!).  We use the
c     straigth insertion method.
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit   none
  
      include            'const.dek'
      include            'sysanal.h'

      integer            i,iorder,j,jidx,k,idx(mxlin),nf(mxlin),n
      logical            sf(mxlin)
      double precision   a(mxlin),fb(mxlin,mxlin),fe(mxlin,mxlin)

      include            'sysanal.com'


c     store the master region array for index sorting

      DO 01 i=1,nlines
        a(i) = f_beg(i,1,1)
 01   CONTINUE  

c     obtain the index order in IDX

      n = nlines
      CALL indexx(n,a,idx)

c     now, loop over each ion transition

      DO 03 iorder=1,norders

c     save the arrays in temporary storage

       DO 05 j=1,nlines
        sf(j) = sf_flag(iorder,j)
        nf(j) = nfind(iorder,j)
        DO 07 k=1,nf(j)
         fb(j,k) = f_beg(j,iorder,k)
         fe(j,k) = f_end(j,iorder,k)
 07     CONTINUE
 05    CONTINUE


c     now move to the IDX location

       DO 15 j=1,nlines
        jidx = idx(j)  
        sf_flag(iorder,j) = sf(jidx)
        nfind(iorder,j)   = nf(jidx) 
        DO 17 k=1,nf(jidx)
         f_beg(j,iorder,k) = fb(jidx,k) 
         f_end(j,iorder,k) = fe(jidx,k)
 17     CONTINUE
 15    CONTINUE

 03   CONTINUE

c     return

      RETURN
      END


c.........................................................................
c

      SUBROUTINE parzeline(nline,ndum,dum,stringline)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      logical           error
      integer           ndum,nline,i,j,ii,istart,iend
      double precision  dum(ndum),value2
      character*80      tmpstr
      character*200     stringline


      DO 11 i=1,ndum
        dum(i) = 0.0d0
 11   CONTINUE

      j     = 0
      i     = 1
      nline = 0


 01    IF (stringline(i:i).ne.' ') then 

         istart = i
         error  = .false.

         DO 15 ii=istart,200
          IF (stringline(ii:ii).ne.' ') then
            iend = ii
          ELSE
            GOTO 02
          END IF
 15      CONTINUE

         nline = j
         RETURN

 02      tmpstr = stringline(istart:iend)
         j      = j + 1
         IF (j.gt.ndum) STOP '(parzeline): J>NDUM'
         dum(j) = value2(tmpstr,error)
         IF (error) dum(j) = 0.0d0
         i = iend + 1
         GOTO 01

       ELSE

         i = i + 1
         GOTO 01
 
       END IF
       

      END

c.........................................................................
c

      SUBROUTINE parzstring(inword,outword1,outword2,error)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      logical             error
      integer             i,k1beg,k1end,k2beg,k2end,lend      
      character*(*)       inword,outword1,outword2


c     find the beginning and ending of each word assuming forms 's1 s2',
c     's1=s2', or 's1 = s2' and allow for leading blanks

      lend  = len(inword)
      k1beg = 0
      k2beg = 0
      k1end = 0
      k1end = 0


c     begin 1st word

      DO 07 i=1,lend
       IF (inword(i:i).ne.' ') then
         k1beg = i
         GOTO 08
       END IF
 07   CONTINUE

c     end first word

 08   DO 09 i=k1beg,lend
       IF ((inword(i:i).eq.' ').OR.(inword(i:i).eq.'=')) then
         k1end = i-1
         GOTO 10
       END IF
 09   CONTINUE

c     begin 2nd word

 10   DO 11 i=k1end+1,lend   
       IF ((inword(i:i).ne.' ').AND.(inword(i:i).ne.'=')) then
         k2beg = i
         GOTO 12
       END IF
 11   CONTINUE

c     end 2nd word

 12   DO 13 i=k2beg,lend
       IF (inword(i:i).eq.' ') then
         k2end = i - 1
         GOTO 14
       END IF
 13   CONTINUE

 14   IF (k1beg*k1end.eq.0) then
        error = .true.
        outword1 = ' '
        outword2 = ' '
        RETURN
      ELSE
        outword1 = inword(k1beg:k1end)
      END IF

      IF (k2beg*k2end.eq.0) then
        outword2 = ' '
      ELSE
        outword2 = inword(k2beg:k2end)
      END IF

      RETURN
      END


c.........................................................................
c

      SUBROUTINE fappend(infile,delim,outfile)   


c     general file appending; this routine seaerched the string "infile"
c     for the location of either the end of the string, defined by the
c     first blank " ", or by the first appearence of the period ".";
c     then the string "delim" is appended after a "."  and stored as
c     string "outfile"

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,k,lend      
      character*(*)       infile,delim,outfile

      lend = len(infile)
      k    = 0

      DO 09 i=1,lend
        k = i
        IF ((infile(i:i).eq.'.').OR.(infile(i:i).eq.' ')) GOTO 10
 09   CONTINUE

 10   IF (delim.eq.' ') then
       outfile = infile(1:k-1)
      ELSE
       IF (infile(k:k).eq.'.') outfile = infile(1:k)//delim
       IF (infile(k:k).eq.' ') outfile = infile(1:k-1)//"."//delim
      END IF


      RETURN
      END


c
c.............................................................................
c

      SUBROUTINE substr1(infile,ionlabel)

c     gets the ion label from the input file, assumes that the
c     information is sandwiched between "." 1 and "." 2

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,k1,k2,lend
      character*(*)       infile,ionlabel


      lend = len(infile)

      k1 = 0
      DO 09 i=1,lend
        k1 = i
        IF (infile(i:i).eq.'.') GOTO 10
 09   CONTINUE


 10   k2 = 0
      DO 11 i=k1+1,lend
        k2 = i
        IF (infile(i:i).eq.'.') GOTO 12
 11   CONTINUE


 12   ionlabel = infile(k1+1:k2-1)

      RETURN
      END

c
c.............................................................................
c

      SUBROUTINE substr2(infile,losindex)

c     gets the los number from the input file, assumes that the
c     information is sandwiched between "." 2 and "." 3

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,k1,k2,k3,lend
      character*(*)       infile,losindex


      lend = len(infile)

      k1 = 0
      DO 09 i=1,lend
        k1 = i
        IF (infile(i:i).eq.'.') GOTO 10
 09   CONTINUE


 10   k2 = 0
      DO 11 i=k1+1,lend
        k2 = i
        IF (infile(i:i).eq.'.') GOTO 12
 11   CONTINUE

 12   k3 = 0
      DO 13 i=k2+1,lend
        k3 = i
        IF (infile(i:i).eq.'.') GOTO 14
 13   CONTINUE


 14   losindex = infile(k2+1:k3-1)

      RETURN
      END

c
c.............................................................................
c

      SUBROUTINE fappend1(infile,delim,outfile)   


c     similar to routine fappend; in this special case, the delimiter of
c     "infile" is assumed to be ".dat"; this is truncated off the
c     "infile" string and the string "delim" is appended following a
c     ":"

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,k,lend      
      character*(*)       infile,delim,outfile

c     lend returns the allocated string length, not the length of the
c     stored string

      lend = len(infile)
      k    = 0

      DO 09 i=1,lend-3
        k = i
        IF (infile(i:i+3).eq.'.dat') GOTO 10  ! k = position of ":"
 09   CONTINUE

 10   outfile = infile(1:k-1)//'.'//delim

      RETURN
      END




c
c.........................................................................
c

      SUBROUTINE fappend2(infile,midstr,delim,outfile)

c     similar to routine fappend1 in that the delimiter of "infile" is
c     assumed to be ".dat"; this is truncated off the "infile" string
c     and then two strings are appended in the form ":midstrdelim"; for
c     flexibility the "." is not automatically placed between "midstr"
c     and "delim"; the use must include the "." in the strong "delim"
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,j,k,lend
      character*(*)       infile,midstr,delim,outfile

c  lend returns the allocated string length, not the length of the
c  stored string

      lend = len(midstr)
      j    = 0

      DO 07 i=1,lend
        j = i
        IF (midstr(i:i).eq.' ') GOTO 08
 07   CONTINUE

 08   lend = len(infile)
      k    = 0

      DO 09 i=1,lend
        k = i
        IF (infile(i:i+3).eq.'.dat') GOTO 10
 09   CONTINUE

 10   outfile = infile(1:k-1)//'.'//midstr(1:j-1)//'.'//delim

      RETURN
      END

c
c.............................................................................
c

      SUBROUTINE sappend(inword,tie,addword,outword)   

c     this is a general string appending routine; the end of "inword" is
c     located by the first blank " "; then, there are two possible
c     appending results; (1) "tie" is blank, in which case "addword" is
c     directly appended to "inword" with no blanks; or (2) "tie" is some
c     string (could be ":" or "-" or something like that, and "addword"
c     is appended to "inword" with "tie" sandwiched inbetween

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,j,k,lend      
      character*(*)       inword,addword,outword,tie

      lend = len(inword)
      k    = 0
      DO 03 i=1,lend
        k = i
        IF (inword(i:i).eq.' ') GOTO 10
 03   CONTINUE

 10   lend = len(tie)
      j    = 0
      DO 05 i=1,lend
        j = i
        IF (tie(i:i).eq.' ') GOTO 11
 05   CONTINUE

 11   outword = inword(1:k-1)//tie(1:j-1)//addword
     

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

      CALL system('/bin/rm -f Mockspec.runlog.anal')

c     remove any old .sysabs files

      CALL system('/bin/rm -f QS0*.sysabs')

c     remove any old .regabs files

      CALL system('/bin/rm -f QS0*.regabs')

      RETURN
      END


c     eof


