c
c.............................................................................
c

      DOUBLE PRECISION FUNCTION REC_Rate(k,j,T)

c     compute the total radiative recombination coefficent at T for
c     species k with ionization stage j; the total recombination
c     coefficient (cm^3/s) is the sum of the radiative, low temperature
c     dielectronic, and high temperature dielectronic rate coefficients

c     IMPORTANT NOTE:
c     for all recombinations, we use the convention of Verner + Yakovlev
c     (1990,A&SS,165,27) [a critical paper!], for which "In all
c     ionization rates, j denotes an initial ionization state, whereas
c     in all recombination rates, j labels the final state."

c     Consider the reaction:  A(k,j) + e- -> A(k,j-1) + E(photon)
c
c     for this reaction, the recombination rate coefficient is
c     alpha(k,j-1), that is, it is indexed by the final ionization level
c     following recombination

c     for example, to compute the recombination reaction rate (number of
c     recombinations per unit volume per second) of CIII to CII, where
c     k=6 and j=3; we have the reaction
c
c     Reaction: A(k,j) + e-  ->  A(k,j-1)
c     ex:       A(6,3) + e-  ->  A(6,2)
c     which is: CIII   + e-  ->  CII
c
c     the reaction rate R(k,j) is calculated as follows:
c
c     R(k,j-1) = ne * n_A(k,j) * alpha(k,j-1)
c     R(6,2)   = ne * n_A(6,3) * alpha(6,2)
c     R(CII)   = ne * n(CIII)  * alpha(CII)

c     one way to think of this is that the rate coefficient is indexed
c     for the ion towards which recombination proceeds

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c


      implicit none
      include           'rates.h'
      integer           k,j
      double precision  T
      double precision  PL,LDR,HDR
      double precision  alpha_rad,alpha_ldr,alpha_hdr
      include           'rates.com'
      include           'com-modes.com'
      include           'com-recomb.com'
      


      alpha_ldr = 0.0d0
      alpha_hdr = 0.0d0

      alpha_rad = PL(k,j,T)             ! radiative
      IF (doLDR) alpha_ldr = LDR(k,j,T) ! low T dielectronic
      IF (doHDR) alpha_hdr = HDR(k,j,T) ! high T dielectronic

      REC_Rate = alpha_rad + alpha_ldr + alpha_hdr

c     store rad and dielectronic components for printing

      alpha_rec_ph(k,j)  = alpha_rad
      alpha_rec_die(k,j) = alpha_ldr + alpha_hdr
 
c     return

      RETURN

      END


c
c.............................................................................
c

      DOUBLE PRECISION FUNCTION PL(k,j,T)

c     compute the radiative recombination coefficent at T towards
c     species k in ionization stage j

c     S. M. V. Aldrovandi & D. Pequignot (1973,A&A,25,137),
c     J. M. Shull & M. Van Steenberg (1982,ApJS,48,95),
c     M. Arnaud & R. Rothenflug (1985,A&AS,60,425).

c     For hydrogenic ions we use the formula given in Sec3.1 of Arnaud +
c     rothenberg (195), taken from Seaton (1959)

c     For non-hydrogenic ions we use the power-law fit to the radiative
c     recombination:
c
c     alpha = A*(T/10**4K)**(-B)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none     
      include           'rates.h'
      integer           j,k,nek
      double precision  T,T4,A,B,Z,lam
      include           'com-recomb.com'



c     initialize to the default null value

      PL = 0.0d0

c     RECOMBINATION TOWARDS HYDROGENIC IONS (e- recombines with fully
c     stripped nucleus; not included in the read in tables)

      nek = k - j + 1

      IF (nek.eq.1) then
       Z   = float(k)
       lam = (1.57890d5/T)*(Z*Z)
       PL  = 5.197d-14*Z*sqrt(lam) * (0.4288d0 + 0.5d0*log(lam) +
     &       0.469*(lam**(-1.0d0/3.0d0)))
       IF (PL.lt.0.0d0) PL = 0.0d0
       RETURN
      END IF


c     RECOMBINATION TOWARDS NON-HYDROGENIC IONS

c     store the fitting parameters for code reading ease; if the k,j
c     combination is off the parameter table then A and B will be null;
c     if this is the case then we use the default null value

      T4 = 1.0d-4*T
      A  = A_pl(k,j)
      B  = B_pl(k,j)

c     the recombination coefficient; check that this is a valid ion

      IF (A.eq.0.0) RETURN

c     if the result is negative , then we null it; this should not
c     happen in this routine, but can happen in the dielectronic
c     routines

      PL = A*(T4**(-B)) 
      IF (PL.lt.0.0d0) PL = 0.0d0

c     return

      RETURN

      END


c
c.............................................................................
c
      DOUBLE PRECISION FUNCTION LDR(k,j,T)

c     compute the low T dielectronic recombination coefficient at T for
c     recombination towards species k in ionization stage j

c     H. Nussbaumer & P. J. Storey (1986,A&A 126, 75)
c     H. Nussbaumer & P. J. Storey (1978,A&AS 64,545)

c     Fit to the of low-temperature dielectronic recombination:

c     alpha = 10**(-12)*(a/t+b+c*t+d*t**2)*t**(-3/2)*exp(-f/t) cm**3/s

c     where t = T[K]/10**4

c     Note that there are two entries for k=8, j=1 (N=8) and 
c                                     for k=8, j=4 (N=5)
c
c     1st -- T < 2*10**4 K
c     2nd -- T > 2*10**4 K

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none     
      include           'rates.h'
      integer           j,k
      double precision  T,T4
      double precision  Tlo,Thi,T4break
      double precision  a,b,c,d,f,Calph
      include           'com-recomb.com'



c     initialize to the default null value

       LDR = 0.0d0

c     the low temperature dielectronic parameter table is quite sparce;
c     thus, we check that we have fitting parameters for this ion; if
c     not, then null and return; we store the number of ionization
c     stages for species k in the array nNstg(k), if j>Nstg(k), then we
c     are off the grid; if the species k is not represented in the
c     parameter table, then nNstg(k)=0

      IF (nNstg(k).eq.0) RETURN
      IF ((j.lt.n1stg(k)).OR.(j.gt.nNstg(k))) RETURN

c     for low temperature dielectronic recombination, the fitting
c     parameters are valid only over the temperature range 1.0e+3 K to
c     6.0e+4 K; if we are outside this range then null the coefficient
c     and return

      Tlo     = 1.0d2 ! test show this is OK, needed for smooth total
      Thi     = 6.0d4

      IF ((T.lt.Tlo).OR.(T.gt.Thi)) RETURN

c     set the normalization constant and normalized temperature

      Calph = 1.0d-12
      T4    = 1.0d-4*T

c     for two of the ions in the table the fitting parameters have a
c     temperature range the demarcation of which is T=2.0e4 K or T4=2;
c     the fitting parameters are dimensioned such that if there are 2
c     temperature ranges for the ion, then we chose the appropriate
c     parameters for the fitting function

      T4break = 2.0d0

c     store the fitting parameters for code reading ease

c     if valid over the full temperature range, then nTreg(j,k)=1

      IF (nTreg(k,j).eq.1) then 
       a = a_ldr(k,j,1)
       b = b_ldr(k,j,1)
       c = c_ldr(k,j,1)
       d = d_ldr(k,j,1)
       f = f_ldr(k,j,1)
      END IF

c     if validity is split over the two ranges demarcated by T4break,
c     then nTreg(j,k)=2; choose and store the correct parameters based
c     upon the temperature

      IF (nTreg(k,j).eq.2) then
       IF (T4.le.T4break) then ! T4bread=2 is logT=4.3
        a = a_ldr(k,j,1)
        b = b_ldr(k,j,1)
        c = c_ldr(k,j,1)
        d = d_ldr(k,j,1)
        f = f_ldr(k,j,1)
       ELSE
        a = a_ldr(k,j,2)
        b = b_ldr(k,j,2)
        c = c_ldr(k,j,2)
        d = d_ldr(k,j,2)
        f = f_ldr(k,j,2)
       END IF
      END IF

c     special low T cut-off check (my doing)

      IF ((f.lt.0.0).AND.(abs(f/T4).gt.0.5d0)) RETURN 

c     the recombination coefficient

      LDR = Calph *(a/T4 + b + c*T4 + d*(T4**2)) * 
     &      T4**(-3.0d0/2.0d0) * exp(-f/T4)

c     even in the range of valid temperatures, there are certain k,j for
c     which the recombination coefficient goes negative (which means it
c     is null); check and null

      IF (LDR.lt.0.0d0) LDR = 0.0d0

c     return

      RETURN

      END


c
c.............................................................................
c
      DOUBLE PRECISION FUNCTION HDR(k,j,T)

c     compute the high temperature dielectronic recombination
c     coefficients for recombination towards species k in ionization
c     stage j

c     S. M. V. Aldrovandi & D. Pequignot (1973, A&A, 25, 137),
c     J. M. Shull & M. Van Steenberg (1982, ApJS, 48, 95),
c     M. Arnaud & R. Rothenflug (1985, A&AS, 60, 425).

c     Fit to the high-temperature dielectronic recombination: Eq.3

c     alpha = A*T**(-3/2)*exp(-T0/T)*[1+B*exp(-T1/T)] cm**3/s

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none     
      include           'rates.h'
      integer           j,k
      double precision  T,Tlo,A,B,T0,T1,Bterm
      include           'com-recomb.com'


c     initilaize to the default null value

      HDR = 0.0d0

c     this is the lower temperature limit of validity; check

      Tlo = 2.0d3

      IF (T.lt.Tlo) RETURN

c     store the fitting parameters for code reading ease; if the k,j
c     combination is off the parameter table then A and B will be null;
c     if this is the case then null the coefficient

      A  = A_hdr(k,j)
      B  = B_hdr(k,j)
      T0 = T0_hdr(k,j)
      T1 = T1_hdr(k,j)

c     the recombination coefficient; check that this is a valid ion

      IF (A.ne.0.0d0) then
       Bterm = 1.0d0+B*exp(-T1/T)
       HDR = A*(T**(-1.50d0))*exp(-T0/T)*Bterm 
      END IF

c     even in the range of valid temperatures, there are certain k,j for
c     which the recombination coefficient goes negative (which means it
c     is null); check and null

      IF (HDR.lt.0.0d0) HDR = 0.0d0

c     return

      RETURN

      END


