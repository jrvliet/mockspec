c
c.........................................................................
c

      SUBROUTINE readCdiA85tab

c     read in the modified table for the fittinf parameters for direct
c     collisional ionization from the work of Arnaud + Rothenflung
c     (1985,A&AS,60,425)

c     the ionization portentials are in units eV; the fitting
c     coefficients, A, B, C, and D are in units 10^-14 cm^2 eV^2, this
c     means we must multiply their Eq 1 by a normalization constant of
c     10^-14 to get the results in cm^2 and by 10^-14/10^-18 = 10^4 to
c     get the result in Mb (10^-18 cm^2)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit none

      include           'rates.h'
      integer           Nrec
      parameter         (Nrec=392)
      integer           i,j,k,s,kold,neold
      integer           z,nel,nsh,lsh
      double precision  ip,aa,bb,cc,dd
      character*80      header
      character*255     tabfile
      include           'rates.com'
      include           'com-coll.com'

c     construct file name

      tabfile = tabpath
      CALL sappend(tabfile,'tab-coll-DI.dat',tabfile)

c     read file

      OPEN(unit=2,file=tabfile,status='old')
      READ(2,*) header

c     hydrogen and helium are single shell; read directly
c     i=1 hydrogen  k=1,j=1,s=1
c     i=2 helium    k=1 j=1 s=1
c     i=3 helium    k=2 j=2 s=1

      s = 1
      DO 11 i=1,3
       READ(2,*,END=32) z,nel,nsh,lsh,ip,aa,bb,cc,dd
       k = z
       j = k - nel + 1
       CALL store(k,j,s,ip,aa,bb,cc,dd)
 11   CONTINUE

c     now cull the file and read in the metals

      READ(2,*,END=32) z,nel,nsh,lsh,ip,aa,bb,cc,dd

 01   k = z
      j = k - nel + 1
      s = 1
      CALL store(k,j,s,ip,aa,bb,cc,dd)

      IF (nel.gt.2) then        ! for nel=1,2 shmax(k,j)=1; don't loop

       neold = nel              ! for change of ion check
       DO s=2,Smax+1            ! loop over possible shells
        READ(2,*,END=32) z,nel,nsh,lsh,ip,aa,bb,cc,dd
        IF (nel.eq.neold) then  ! same ion; k,j don't change
         CALL store(k,j,s,ip,aa,bb,cc,dd)
        ELSE                    ! new ion; k a/o j change; s=1
         GOTO 01 
        END IF
       END DO

      ELSE                      ! single shell case; k a/o j change; s=1 

       READ(2,*,END=32) z,nel,nsh,lsh,ip,aa,bb,cc,dd
       GOTO 01

      END IF


 32   CLOSE(unit=2)


c     DEBUG?

C      DO 41 k=1,Imax
C       DO 42 j=1,k
C        DO 43 s=1,shmax(k,j)
C         WRITE(99,600) k,j,s,IE(k,j,s),A(k,j,s),
C     &                 B(k,j,s),C(k,j,s),D(k,j,s)  
C 43     CONTINUE
C 42    CONTINUE
C 41   CONTINUE
C 600  FORMAT(1x,3i3,5(1x,f9.3))

      RETURN

      END


c
c.........................................................................
c

      SUBROUTINE store(k,j,s,ip,aa,bb,cc,dd)

c     this routine simply stores the read in line into the common block
c     array and matrices for the kiven k,j,s combination

c     IE(k,j,s)     = threshold energy for shell s for k,j
c     A(k,j,s)      = Arnaud A coefficent for k,j,s
c     B(k,j,s)      = Arnaud B coefficent for k,j,s
c     C(k,j,s)      = Arnaud C coefficent for k,j,s
c     D(k,j,s)      = Arnaud D coefficent for k,j,s

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           j,k,s
      double precision  ip,aa,bb,cc,dd

      include           'com-coll.com'


c     store the fitting parameters and shell data

      IE(k,j,s) = ip
      A(k,j,s)  = aa
      B(k,j,s)  = bb
      C(k,j,s)  = cc
      D(k,j,s)  = dd

c     increment the shell counter for k,j

      shmax(k,j) = shmax(k,j) + 1

      RETURN
      END

