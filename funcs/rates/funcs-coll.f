c
c
c.............................................................................
c

      DOUBLE PRECISION FUNCTION CDI_Rate(k,j,T)

c     compute the total direct collisional rate coefficient for ion k,j
c     at temperature T

c     this function evaluates the 2nd equation (it is unnumbered) on
c     page 426 from Arnaud + Rothenflung (1985,A&AS,60,425)

c     each term of the sum is the partial direct collisional rate
c     coefficient for ion k,j for shell=s; we perform the sum here, so
c     it is the total direct collisional rate coefficient for ion k,j
c     that is returned

c     Arnaud + Rothenflung index the shells with index j; we index the
c     shells with s (where k is the species atomic number and j is the
c     ionization stage)

c     the units of the returned direct collisional rate coefficient is
c     cm^3/s

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit none

      include           'rates.h'
      integer           k,j,s
      double precision  F,x
      double precision  T,kT,sum
      double precision  AA,BB,CC,DD
      include           'com-coll.com'


c     set to the default value

      CDI_Rate = 0.0d0

c     convert the temperature to energy

      kT = 8.61738573d-5*T   ! eV

c     perform the sum over the shells (s) for k,j; the function F(x_j)
c     [as denoted by Arnaud + Rothenflung] is the function call
c     F(x,AA,BB,CC,DD) in the loop; the conidtion for x<80 avoids
c     underflow with the exponential function

      sum = 0.0d0
      DO 11 s=1,shmax(k,j)
       x   = IE(k,j,s)/kT   ! x  denoted x_j by Arnaud + Rothenflung 
       AA  = A(k,j,s)       ! AA denoted A_j 
       BB  = B(k,j,s)       ! BB denoted B_j 
       CC  = C(k,j,s)       ! CC denoted C_j 
       DD  = D(k,j,s)       ! DD denoted D_j 
       IF (x.ne.0.0.AND.x.lt.80.0d0) then
        sum = sum + F(x,AA,BB,CC,DD)*exp(-x)/x
      END IF
 11   CONTINUE

      CDI_Rate = 6.69d-7*sum/(kT**1.5d0)

      IF (CDI_Rate.lt.0.0) CDI_Rate = 0.0d0

      RETURN

      END


c
c.............................................................................
c

      DOUBLE PRECISION FUNCTION CEA_Rate(k,j,T)

c     compute the excitation-autoionization collisional rate coefficient
c     for ion k,j at temperature T

c     this function evaluates the C_EA equations (they are unnumbered)
c     from Appendix A (page 435) from Arnaud + Rothenflung
c     (1985,A&AS,60,425); the particular equation depends upon the
c     number of bound electrons in the ion (the isoelectronic series)

c     the units of the returned excitation-autoionization rate
c     coefficient is cm^3/s

c     I was able to organize the equations so that they all took the
c     same form; but this means that the normalization constants do not
c     appear exactly as written in the paper; note that Gy=e^-y*G(y)
c     where G(y) is as given in the paper

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit none

      include           'rates.h'
      integer           k,j,nel
      double precision  f1
      double precision  Iea,y,a0,x,zeff,f1y,Gy
      double precision  T,kT,Z
      double precision  b0
      include           'rates.com'
      include           'com-coll.com'


c     set to the default null value

      CEA_Rate = 0.0d0

c     obtain the number of electrons (which defines the series); if this
c     species has no data in the table, which means nek(k,j)=0, then
c     bail

       nel = k - j + 1

c      IF (nel.eq.0) RETURN

c     the first three species have no excitation-autoionization physics;
c     except for four special cases, ions with 19 electrons have no AE
c     process

      IF (k.le.3) RETURN
      IF (nel.eq.3) GOTO 01  ! lithium series good for all k

      IF (nel.gt.18) then
       IF (k.eq.20.AND.nel.eq.20) GOTO 01 ! Ca^0  (CaI)  (k,j,nel)=(20,1,20)
       IF (k.eq.20.AND.nel.eq.19) GOTO 01 ! Ca^+1 (CaII) (k,j,nel)=(20,2,19)
       IF (k.eq.26.AND.nel.eq.23) GOTO 01 ! Fe^+3 (FeIV) (k,j,nel)=(26,4,23) 
       IF (k.eq.26.AND.nel.eq.22) GOTO 01 ! Fe^+4 (FeV)  (k,j,nel)=(26,5,22)
       RETURN
      END IF

c     there are a several sequences for which the excitation-autoionization
c     is null or negligible

      IF (nel.eq.1)  RETURN ! hydrogen  sequence (1s^1)
      IF (nel.eq.2)  RETURN ! helium    sequence (1s^2)
      IF (nel.eq.4)  RETURN ! beryllium sequence (2s^2)
      IF (nel.eq.5)  RETURN ! boron     sequence (2p^1)
      IF (nel.eq.6)  RETURN ! carbon    sequence (2p^2)
      IF (nel.eq.7)  RETURN ! nitrogen  sequence (2p^3)
      IF (nel.eq.8)  RETURN ! oxygen    sequence (2p^4)
      IF (nel.eq.9)  RETURN ! florine   sequence (2p^5)
      IF (nel.eq.10) RETURN ! neon      sequence (2p^6)
      IF (nel.eq.17) RETURN ! chlorine  sequence (3p^5)
      IF (nel.eq.18) RETURN ! argon     sequence (3p^6)

c     apart from lithium, all series for k>28 DNE

      IF (k.gt.28) RETURN    

c     neutral sodium and neutral magnesium have no
c     excitation-autoionization process

      IF (k.eq.11.AND.nel.eq.11) RETURN  ! neutral sodium
      IF (k.eq.12.AND.nel.eq.12) RETURN  ! neutral magnesium

c     between aluminum (k=13) and sulfur (k=16), the aluminum,
c     magnesium, silicon, phosphorous, and sulfur sequences are
c     negligible; but lithium is OK

      IF (k.ge.13.AND.k.le.16.AND.nel.ge.12) RETURN

c     for chlorine (k=17), the aluminum, magnesium, silicon,
c     phosphorous, and sulfur sequences are negligible; but lithium is
c     OK

      IF (k.eq.17.AND.nel.le.16.AND.nel.ge.11) RETURN

c
c     BELOW SEQUENCES ARE GOOD TO GO
c

c     convert the temperature to energy and obtain the species nuclear
c     charge Z

 01   kT  = 8.61738573d-5*T     ! eV
      Z   = float(k)

c     based upon nel, set up the parameters for the computation of the
c     rate coefficient

c     lithium sequence (1p^1); nel=3 (Z>=4); this sequence includes the
c     famous resonant doublet transition ions BeII, CIV, NV, OVI, and
c     NeVIII and other not so well observed very high ions; the
c     normalizing correction term b0=1 applies for only the 1s-2p
c     transition; non-unity b0 accounts as an approximation for the
c     other transitions and has an ion dependence; the fiducual b0=1.2

      IF ((nel.eq.3).AND.(k.ge.4)) then 
       Iea  = 13.6d0*((Z-0.835d0)**2 - 0.25d0*(Z-1.62d0)**2)
       y    = Iea/kT
       IF (y.gt.80.0) RETURN  ! avoid underflow of exponential
       b0   = 1.20d0
       IF (k.eq.6) b0 = 0.6d0   ! carbon   (CIV)    
       IF (k.eq.7) b0 = 0.8d0   ! nitrogen (NV)
       IF (k.eq.8) b0 = 1.25d0  ! oxygen   (OVI)
       a0   = (1.60d-7)*b0
       x    = 1.0d0 + (2.0d-4)*Z*Z*Z   ! branching ratio
       zeff = (Z-0.43d0)    
       f1y  = f1(y)
       Gy   = exp(-y)*(2.22d0*f1y + 0.67d0*(1.0d0-y*f1y) + 
     &        0.49d0*y*f1y +1.2d0*y*(1.0d0-y*f1y)) 
       CEA_Rate = a0*Gy/(sqrt(kT)*x*(zeff**2))
       RETURN
      END IF

c     low Z sodium sequence (3s^1); nel=11 (12<=Z<=16) magnesium to
c     sulfer ions; this low Z sequence includes the famous resonant
c     doublet transition ions MgII, AlIII, SiIV, PV,and SVI
 
       IF ((nel.eq.11).AND.(k.ge.12).AND.(k.le.16)) then 
       Iea  = 26.0d0*(Z-10.0d0)
       y    = Iea/kT
       IF (y.gt.80.0) RETURN  ! avoid underflow of exponential
       a0   = 1.8732d-9
       x    = 1.0d0/Iea
       zeff = (Z-11.0d0)**0.35d0
       f1y  = f1(y)
       Gy   = exp(-y)*(1.0d0-y*f1y)
       CEA_Rate = a0*Gy/(sqrt(kT)*x*(zeff**2))
       RETURN
      END IF

c     high Z sodium sequence (3s^1); nel=11 (18<=Z<=28) argon to nickel
c     ions; this high Z sequence includes the not as famous resonant
c     doublet transition ions ArVIII, CaX, TiXII, FeXVI, etc

      IF ((nel.eq.11).AND.(k.ge.18).AND.(k.le.28)) then 
       Iea  = 11.0d0*(Z-10.0d0)**1.5d0
       y    = Iea/kT
       IF (y.gt.80.0) RETURN  ! avoid underflow of exponential
       a0   = 8.697d-7
       x    = 1.0d0/Iea
       zeff = (Z-10.0d0)**1.865d0
       f1y  = f1(y)
       Gy   = exp(-y)*(1.0d0-0.5d0*y*(1.0d0-y*(1.0d0-y*f1y)))
       CEA_Rate = a0*Gy/(sqrt(kT)*x*(zeff**2))
       RETURN
      END IF

c     high Z (18<=Z<=28) magnesium (3s^2) through sulfer 3p^4 sequences
c     Mg nel=12 ArVII,KVIII,CaIX,ScX,TiXI,VXII,CrXIII,MnXIV,FeXV.CoXVI,NiXVII
c     Al nel=13 ArVI,KVII,CaVIII,ScIX,TiX,VXI,CrXII,MnXIII,FeXIV.CoXV,NiXVI
c     Si nel=14 ArV,KVI,CaVII,ScVIII,TiIX,VX,CrXI,MnXII,FeXIII.CoXIV,NiXV
c     P  nel=15 ArIV,KV,CaVI,ScVII,TiVIII,VIX,CrX,MnXI,FeXII.CoXIII,NiXIV
c     S  nel=16 ArV,KIV,CaV,ScVI,TiVII,VVIII,CrIX,MnX,FeXI.CoXII,NiXIII

      IF ((nel.ge.12).AND.(nel.le.16).AND.(k.ge.18).AND.(k.le.28)) then 
       IF (nel.eq.12) Iea = 10.3d0*(Z-10.0d0)**1.52d0
       IF (nel.eq.13) Iea = 18.0d0*(Z-11.0d0)**1.33d0
       IF (nel.eq.14) Iea = 18.4d0*(Z-12.0d0)**1.36d0
       IF (nel.eq.15) Iea = 23.7d0*(Z-13.0d0)**1.29d0
       IF (nel.eq.16) Iea = 40.0d0*(Z-14.0d0)**1.10d0
       y    = Iea/kT
       IF (y.gt.80.0) RETURN  ! avoid underflow of exponential
       a0   = 2.676d-5
       x    = 1.0d0
       zeff = Z
       f1y  = f1(y)
       Gy   = exp(-y)*(1.0d0-0.5d0*y*(1.0d0-y*(1.0d0-y*f1y)))
       CEA_Rate = a0*Gy/(sqrt(kT)*x*(zeff**2))
       RETURN
      END IF

c     OK, it must be one of the four special cases

      IF ((k.eq.20).AND.(nel.eq.20)) then  ! Ca^0  (CaI)  (k,j,nel)=(20,1,20)
       Iea  = 25.0d0
       y    = Iea/kT
       IF (y.gt.80.0) RETURN  ! avoid underflow of exponential
       a0   = 4.014d-9
       x    = 1.0d0/Iea
       zeff = 1.0d0
       f1y  = f1(y)
       Gy   = exp(-y)*(1.0d0+1.12d0*f1y)
       CEA_Rate = a0*Gy/(sqrt(kT)*x*(zeff**2))
       RETURN
      END IF

      IF ((k.eq.20).AND.(nel.eq.19)) then  ! Ca^+1 (CaII) (k,j,nel)=(20,2,19)
       Iea  = 29.0d0
       y    = Iea/kT
       IF (y.gt.80.0) RETURN  ! avoid underflow of exponential
       a0   = 6.5562d-9
       x    = 1.0d0/Iea
       zeff = 1.0d0
       f1y  = f1(y)
       Gy   = exp(-y)*(1.0d0+1.12d0*f1y)
       CEA_Rate = a0*Gy/(sqrt(kT)*x*(zeff**2))
       RETURN
      END IF

      IF ((k.eq.26).AND.(nel.eq.23)) then  ! Fe^+3 (FeIV) (k,j,nel)=(26,4,23) 
       Iea  = 60.0d0
       y    = Iea/kT
       IF (y.gt.80.0) RETURN  ! avoid underflow of exponential
       a0   = 1.2042d-9
       x    = 1.0d0/Iea
       zeff = 1.0d0
       f1y  = f1(y)
       Gy   = exp(-y)*(1.0d0+f1y)
       CEA_Rate = a0*Gy/(sqrt(kT)*x*(zeff**2))
       RETURN
      END IF

      IF ((k.eq.26).AND.(nel.eq.22)) then  ! Fe^+4 (FeV)  (k,j,nel)=(26,5,22)
       Iea  = 73.0d0
       y    = Iea/kT
       IF (y.gt.80.0) RETURN  ! avoid underflow of exponential
       a0   = 3.345d-9
       x    = 1.0d0/Iea
       zeff = 1.0d0
       f1y  = f1(y)
       Gy   = exp(-y)*(1.0d0+f1y)
       CEA_Rate = a0*Gy/(sqrt(kT)*x*(zeff**2))
       RETURN
      END IF

c     if we get here, then we are in a bad way

      WRITE(6,*) 'FAIL: FUNCTION CEA_Rate- no trap set'
      WRITE(6,*) '(k,j,ne)=',k,j,nel,sequence(nel)(1:30)
      STOP

      END

c
c.............................................................................
c

      DOUBLE PRECISION FUNCTION F(x,a,b,c,d)

c     this function evaluates the function F(x) from Arnaud +
c     Rothenflung (1985,A&AS,60,425) 

c     used by CDI_Rate

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      double precision  f1,f2,f1x,x
      double precision  a,b,c,d

c     evaluate f1(x) only once

      f1x = f1(x)

c     compute F(x)

      F = a*(1.0d0-x*f1x) + b*(1.0d0+x-x*(2.0d0+x)*f1x) + 
     &    c*f1x + d*x*f2(x)

      RETURN

      END

c
c.............................................................................
c

      DOUBLE PRECISION FUNCTION f1(x)

c     this function evaluates the function f1(x) from Arnaud +
c     Rothenflung (1985,A&AS,60,425) as outline in Appendix B

c     this is supposed to be good to better than 1% accuracy for all x

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      double precision  x,a


c     for x<=0.02

      IF (x.le.0.02d0) then
       f1 = (x-0.5772d0-log(x))*exp(x)
       RETURN
      END IF

c     for x>0.02 and x<1.5 (a=-0.5)

      IF ((x.gt.0.02d0).AND.(x.lt.1.5d0)) then
       a  = -0.5d0
       f1 = log((x+1.0d0)/x) -
     &      (0.36d0+0.03d0*(x+0.01d0)**a)/((x+1.0d0)**2)
       RETURN
      END IF

c     for x>=1.5 and x<10 (a=+0.5)

      IF ((x.ge.1.5d0).AND.(x.lt.10.0d0)) then
       a  = 0.5d0
       f1 = log((x+1.0d0)/x) -
     &      (0.36d0+0.03d0*(x+0.01d0)**a)/((x+1.0d0)**2)
       RETURN
      END IF

c     for x>=10

      f1 = ( 1.0d0 - 1.0d0/x + 2.0d0/(x*x) - 6.0d0/(x*x*x) +
     &       24.0d0/(x*x*x*x) ) / x 


      RETURN

      END

c
c.............................................................................
c

      DOUBLE PRECISION FUNCTION f2(x)

c     this function evaluates the function f2(x) from Arnaud +
c     Rothenflung (1985,A&AS,60,425) as outline in Appendix B

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer           i,ji
      double precision  x,Px,Qx
      double precision  p,q
      COMMON /pxqx/     p(13),q(14) 


c     according to Verner (as noted in his ci.txt file) this is the
c     prefered evaluation of f2(x) for x<0.27

      IF (x.lt.0.27d0) then
       f2 = 0.5d0*log(x)*log(x) + 0.577722d0*log(x) + 1.0d0
       RETURN
      END IF

c     for all other x

c     initialize the 0th terms in the sums

      Px = 1.0d0
      Qx = 1.0d0

c     sum the P(x) term

      DO 11 i=1,13
       ji = float(i)
       Px = Px + p(i)*(x**(-ji)) 
 11   CONTINUE

c     sum the Q(x) term

      DO 15 i=1,14
       ji = float(i)
       Qx = Qx + q(i)*(x**(-ji)) 
 15   CONTINUE

      f2 = Px/(x*x*Qx)

      RETURN

      END

c
c.............................................................................
c

      BLOCK DATA f2data

c     this data block stores the p_j and q_j coefficients for the
c     computations of the sums P(x) and Q(x) for the evaluation of f2(x)
c     for x>=0.27; the coefficeints are given in Appendix B of Arnaud +
c     Rothenflung (1985,A&AS,60,425); since the first terms of the sums
c     are 1.0, we index from 1-13 and 1-14 (instead of 0-13 and 0-14),
c     and explcitly evaluate the 0th term in the above function call

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit none
      integer           i
      double precision  p,q
      COMMON /pxqx/     p(13),q(14) 

      DATA (p(i),i=1,13) /2.1658D+02,2.0336D+04,1.0911D+06,
     &                    3.7114D+07,8.3963D+08,1.2889D+10,
     &                    1.3449D+11,9.4002D+11,4.2571D+12,
     &                    1.1743D+13,1.7549D+13,1.0806D+13,
     &                    4.9776D+11/

      DATA (q(i),i=1,14) /2.1958D+02,2.0984D+04,1.1517D+06,
     &                    4.0349D+07,9.4900D+08,1.5345D+10,
     &                    1.7182D+11,1.3249D+12,6.9071D+12,
     &                    2.3531D+13,4.9432D+13,5.7760D+13,
     &                    3.0225D+13,3.3641D+12/

      END
