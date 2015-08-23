c     PACKAGE: Mockspec
c     PROGRAM: specsynth
c     MODULE:  strings.f
c
c     
c     DESCRIPTION
c     this module contains string manipulation routines that are
c     specific to specsynth
c
c
c     this file contains:
c     SUBROUTINE sappend
c     SUBROUTINE fappend1
c     SUBROUTINE fappend2
c
c
c.............................................................................
c
c

      SUBROUTINE sappend(instring,appstring,outstring)   

c     append appstring on instring and output the result as outstring

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
 
      integer             i,k,lend      
      character*(*)       instring,appstring,outstring


c     lend returns the allocated string length, not the length of the
c     stored string

      lend = len(instring)
      k    = 0
      do 09 i=1,lend
       k = i
       if (instring(i:i).eq.' ') goto 10
 09   continue

 10   outstring = instring(1:k-1)//appstring      

c     normal exit

      return

      end


c
c.............................................................................
c

      SUBROUTINE fappend1(infile,delim,outfile)   

c     file appending subroutine, assumes we are replacing ".dat"

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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


c
c.............................................................................
c

      SUBROUTINE fappend2(infile,midstr,delim,outfile)   

c     file appending subroutine, assumes we are replacing ":*" delim
c     includes "."

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,j,k,lend      
      character*(*)       infile,midstr,delim,outfile

c     lend returns the allocated string length, not the length of the
c     stored string

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

 10   outfile = infile(1:k-1)//'.'//midstr(1:j-1)//delim

      RETURN
      END

c
c.............................................................................
c

      SUBROUTINE substr1(infile,ionlabel)   

c     file appending subroutine, assumes we are replacing ".dat"

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


