c
c.............................................................................
c

      SUBROUTINE readAugerKM93

c     read in the modified table for the fitting parameters for direct
c     collisional ionization from the work of Kaastra + Mewe
c     (1993.A&ASS,97,443)

c     we convert the Xray shells indices of Kaastra to the
c     photoionization shell indices of Verner because we will need to
c     compute the Auger rates from the partial photoionization rates for
c     each shell and these rates are computed using the Verner indices

c     the translation is as follows:
c
c     Kaastra                Verner
c     idx  Xray  Term        s     shell  g
c      1   K     1s_1/2      1     1s     
c      2   L1    2s_1/2      2     2s
c      3   L2    2p^5_1/2    3     2p     1  
c      4   L3    2p_5_3/2    3     2p     2   
c      5   M1    3s_1/2      4     3s
c      6   M2    3p5_1/2     5     3p     1
c      7   M3    3p5_3/2     5     3p     2
c      8   M4    3d9_3/2     6     3d
c      9   M5    3d9_5/2     6     3d
c     10   N1    4s          7     4s
c
c     where g is the statistical weight for the averaging; Kaastra's
c     table goes no higher than idx=7

c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit none

      include           'rates.h'
      integer           Nrec
      parameter         (Nrec=1090)
      integer           i,ii,j,k,n,s,nline,Xshnum,z
      integer           jAugktot(Imax),nAugshtot(Imax,Imax),
     &                  Xshell(Imax,Imax,Nshmax)
      double precision  W(Imax,Imax,Nshmax,10)
      double precision  pr(10),E,EA,eps
      character*80      header
      character*255     tabfile
      include           'rates.com'
      include           'com-auger.com'

c     zero counters and sums [critical]

      DO k=1,Imax
       jAugktot(k) = 0
       DO j=1,Imax
        nAugshtot(k,j) = 0
        DO s=1,Nshmax
         Xshell(k,j,s) = 0
        ENDDO
       ENDDO
      ENDDO

c     construct file name

      tabfile = tabpath
      CALL sappend(tabfile,'tab-auger-Wijk.dat',tabfile)

c     read file

      OPEN(unit=2,file=tabfile,status='old')
      READ(2,*) header
      DO 11 i=1,Nrec
       READ(2,*,END=12) Nline,z,j,Xshnum,E,EA,eps,(pr(ii),ii=1,10)
       k = z
       CALL store_auger(k,j,Xshnum,pr,W,nAugshtot,jAugktot)
 11   CONTINUE
 12   CLOSE(unit=2)

c     translate the Xray shell numbers to Verner's shell numbers so we
c     can calculate the Auger rates using the photoionization cross
c     sections from Verner

c     the translation and weighting is given in the header comments of
c     this routine

      DO 21 k=1,Imax
       DO 22 j=1,jAugktot(k)
        DO 23 n=1,nAugshtot(k,j)

         IF (Xshell(k,j,n).eq.1) then
          s = 1
          DO 31 i=1,10
           W_Auger(k,j,s,i) = W(k,j,n,i) 
 31       CONTINUE
         END IF

         IF (Xshell(k,j,n).eq.2) then
          s = 2
          DO 32 i=1,10
           W_Auger(k,j,s,i) = W(k,j,n,i) 
 32       CONTINUE
         END IF

         IF (Xshell(k,j,n).eq.3) then
          s = 3
          DO 33 i=1,10
           W_Auger(k,j,s,i) = (W(k,j,n,i) + 2.0d0*W(k,j,n+1,i))/3.0d0
 33       CONTINUE
         END IF

         IF (Xshell(k,j,n).eq.5) then
          s = 4
          DO 34 i=1,10
           W_Auger(k,j,s,i) = W(k,j,n,i) 
 34       CONTINUE
         END IF

         IF (Xshell(k,j,n).eq.6) then
          s = 5
          DO 35 i=1,10
           W_Auger(k,j,s,i) = (W(k,j,n,i) + 2.0d0*W(k,j,n+1,i))/3.0d0
 35       CONTINUE
         END IF

 23     CONTINUE
 22    CONTINUE
 21   CONTINUE


      RETURN

      END


c
c.............................................................................
c

      SUBROUTINE store_auger(k,j,Xshnum,pr,W,nAugshtot,jAugktot)

c     this routine simply stores the read in line into the common block
c     array and matrices for the given k,j combination

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           i,j,k,n,Xshnum
      integer           jAugktot(Imax),nAugshtot(Imax,Imax),
     &                  Xshell(Imax,Imax,Nshmax)
      double precision  W(Imax,Imax,Nshmax,10)
      double precision  pr(10)
      include           'com-auger.com'


c     store the fitting parameters and subshell data

      n  = nAugshtot(k,j) + 1

      Xshell(k,j,n) = Xshnum

      DO 11 i=1,10
       W(k,j,n,i) = 1.0d-4*pr(i)
 11   CONTINUE


c     save the shell counters for k,j

      nAugshtot(k,j) = n
      jAugktot(k)    = max(jAugktot(k),j)

      RETURN
      END

