c     PACKAGE: Mockspec
c     PROGRAM: general module for multiple programs in the package
c     MODULE:  parze.f
c
c     this file contains:
c     SUBROUTINE parzehdr
c     SUBROUTINE parzeline
c     DBL FUNCTION value2
c
c
c.............................................................................
c

      SUBROUTINE parzehdr(ndum,dum,stringline)

c     general string parzing subroutine for getting header quantities of
c     files; the single input stringline is parzed into ndum individual
c     strings stored in array 1D dum; the individual strings are parzed
c     under the assumption that they are separated by a blank character
c     in the string stringline

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
 
      implicit none

      logical           error
      integer           ndum,nline,i,j,ii,istart,iend
      character*80      tmpstr,dum(ndum)
      character*250     stringline


      DO 11 i=1,ndum
        dum(i) = ''
 11   CONTINUE

      j     = 0
      i     = 1
      nline = 0


 01       IF (stringline(i:i).ne.' ') then 

         istart = i
         error  = .false.

         DO 15 ii=istart,250
          IF (stringline(ii:ii).ne.' ') then
            iend = ii
          ELSE
            GOTO 02
          END IF
 15     CONTINUE

         nline = j
         RETURN

 02      tmpstr = stringline(istart:iend)
         j      = j + 1
         IF (j.gt.ndum) THEN
           WRITE(6,*) 'ERROR(parzehdr): too many fields in stringline'
           WRITE(6,*) ' ',j,' strings parzed but only ',ndum,' were'
           WRITE(6,*) ' expected based upon passed variable NDUM'
           STOP 
         END IF
         dum(j) = stringline(istart:iend)
         i = iend + 1
         GOTO 01

       ELSE

         i = i + 1
         GOTO 01
 
       END IF

      END


c
c.............................................................................
c

      SUBROUTINE parzeline(ndum,dum,stringline)

c     identical to parzehdr except that the fields of strongline are
c     assume to be numeric and are converted to reals using routine
c     value2; if numerical conversion fails then we null the value
      
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      logical           error
      integer           ndum,nline,i,j,ii,istart,iend
      double precision  dum(ndum),value2
      character*80      tmpstr
      character*250     stringline


      DO 11 i=1,ndum
        dum(i) = 0.0d0
 11   CONTINUE

      j     = 0
      i     = 1
      nline = 0


 01       IF (stringline(i:i).ne.' ') then 

         istart = i
         error  = .false.

         DO 15 ii=istart,250
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
         IF (j.gt.ndum) THEN
           WRITE(6,*) 'ERROR(parzeline): too many fields in stringline'
           WRITE(6,*) ' ',j,' strings parzed but only ',ndum,' were'
           WRITE(6,*) ' expected based upon passed variable NDUM'
           STOP 
         END IF
         dum(j) = value2(tmpstr,error)
         IF (error) dum(j) = 0.0d0
         i = iend + 1
         GOTO 01

       ELSE

         i = i + 1
         GOTO 01
 
       END IF

      END


c
c.............................................................................
c

      double precision function value2(string,err)

c     function value2 converts strings to numerical values if it can
c     when it fails, it returns ERR=.true.  Code obtained unmodified
c     from Frank Timmes (fxt)


c     this routine takes the character string number and converts it to
c     a real. if trouble is encountered during the conversion, this
c     routine returns the logical err as .true.

c     required routines: none

c     date code: 28may90 fxt

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      logical          pflag,err
      character*(*)    string
      character*1      plus,minus,decmal,blank,se,sd,se1,sd1
      integer          noblnk,long,ipoint,power,psign,iten,j,z,i
      double precision x,sign,factor,rten,temp
      parameter        (plus = '+'  , minus = '-' , decmal = '.'   ,
     1                  blank = ' ' , se = 'e'    , sd = 'd'       ,
     2                  se1 = 'E'   , sd1 = 'D'   , rten =  10.0,
     3                  iten = 10                                   )
c..
c..initialize
      err    =  .false.
      x      =  0.0d0
      sign   =  1.0d0
      factor =  rten
      pflag  =  .false.
      noblnk =  0
      power  =  0
      psign  =  1
      long   =  len(string)
c..
c..remove leading blanks and the sign
      do 10 z = 1,7
       noblnk = noblnk + 1
       if ( string(noblnk:noblnk) .eq. blank) then
        if (noblnk .gt. 6 ) goto 1000
       else
        if (string(noblnk:noblnk) .eq. plus) then
         noblnk = noblnk + 1
        else if (string(noblnk:noblnk) .eq. minus) then
         noblnk = noblnk + 1
         sign =  -1.0d0
        end if
        go to 100
       end if
 10   continue
c..
c..main number conversion loop
 100  do 200 i = noblnk,long
       ipoint = i + 1
c..
c..if blank character then we are done
       if ( string(i:i) .eq. blank ) then
        x = x * sign
        value2 = x 
        return


c..
c..if it is an exponent process it
       else if (string(i:i).eq.se  .or. string(i:i).eq.sd .or.
     1          string(i:i).eq.se1 .or. string(i:i).eq.sd1   ) then
        if (x .eq. 0.0d0 .and. ipoint.eq.2)     x = 1.0d0
        if (sign .eq. -1.0d0 .and. ipoint.eq.3) x = 1.0d0
        if (string(ipoint:ipoint) .eq. plus) ipoint = ipoint + 1
        if (string(ipoint:ipoint) .eq. minus) then
         ipoint = ipoint + 1
         psign = -1
        end if
        do 150 z = ipoint,long
         if (string(z:z) .eq. blank)  then
          x = sign * x * rten**(power*psign)
          value2 = x
          return
         else
          j = ichar(string(z:z)) - 48
          if ( (j.lt.0) .or. (j.gt.9) ) go to 1000
          power= (power * iten)  + j
         end if
 150    continue
c..
c..if it is a number process it
       else if (string(i:i) .ne. decmal) then
        j = ichar(string(i:i)) - 48
        if ( (j.lt.0) .or. (j.gt.9) ) go to 1000
        if (.not.(pflag) ) then
         x = (x*rten) + j
        else
         temp   = j
         x      = x + (temp/factor)
         factor = factor * rten
         go to 200
        end if
c..
c..must be a decimal point
       else
        if (pflag) go to 1000
        pflag = .true.
       end if
 200  continue
c..
c..errors
 1000 err = .true.
      return
      end

