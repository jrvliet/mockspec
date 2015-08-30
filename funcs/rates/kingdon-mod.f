*     Code by Jim Kingdon, in collaboration with G.J. Ferland
c     also some code by Chris Churchill

c     modified by Chris Churchill for double precision, additional
c     comments (use c; original comments by Kingdon use *) and for
c     passing the temperature in the argument list instead of via a
c     common block; also several other mods for clarity of reading


******************************************************************************

      DOUBLE PRECISION FUNCTION HCTRecomb(nelem,ion,Te)

c     returns the hydrogen charge exchange recombination rate
c     coefficient (cm^3/s) using the work of Kingdon + Ferland
c     (1996,ApJS,106,205)

c     HCTRecomb is returned in units cm^3/s

*     ion is stage of ionization, 2 for the ion going to the atom
*     nelem is atomic number of element, 2 up to 30
*     Example:  O+ + H => O + H+ is HCTRecomb(8,2,Te)

      implicit none
      integer           ipIon,ion,nelem
      double precision  Te,Tmin,Tmax,T4 
      double precision  a,b,c,d

      double precision  H1r
      common/CTHRecomb/ H1r(6,4,30)



      HCTRecomb = 0.0d0

c     if off the element charge return; if this is the neutral ion, then
c     no CE, return

      IF (nelem.lt.2.OR.nelem.gt.30) RETURN
      IF (ion.eq.1) RETURN

c     note the translation to ion-1 !!!

      ipIon = ion - 1

c     for ions with j>4 we use the statistical results of Ferland,
c     Korista, & Verner (1997,ApJ,481,L115) Eq. 3; the constant 1.92E-9
c     is the reduced mass of the ion in terms of the electron mass for
c     hydrogen interactions; this is the limiting value for large mass
c     ions (i.e., iron); NOTE: they state that a factor of X uncertainty
c     in this constant yields a factor of uncertainty of X^n in the
c     resulting ionization fractions for Ipion=n+5 (n>=1)

      IF (ipIon.gt.4) then  ! use statistical charge transfer for ion > 4
        HCTRecomb = 1.92e-9 * float(ipIon)
        RETURN
      END IF

c     if the first fitting parameter is null, then no data for this ion;
c     hoping that this is because the rate is negligible!

      IF (H1r(1,ipIon,nelem).eq.0.0) RETURN

c     to stay within the allowed T range, this sets the T to be used in
c     the fitting function equal to the limiting Tmin or limiting Tmax;
c     this avoids setting the rate coefficient to zero while avoiding
c     extrapolation

      Tmin  = H1r(5,ipIon,nelem)
      Tmax  = H1r(6,ipIon,nelem)

      IF (Te.lt.Tmin) RETURN
      IF (Te.gt.Tmax) then
       T4 = Tmax * 1.0d-4
      ELSE
       T4 = Te * 1.0d-4
      END IF


c     compute the rate from the fitting function for this T

      a = H1r(1,ipIon,nelem)
      b = H1r(2,ipIon,nelem)
      c = H1r(3,ipIon,nelem)
      d = H1r(4,ipIon,nelem)

      HCTRecomb = 1.0d-9 * a*(T4**b)*(1.0d0 + c*exp(d*T4))

      IF (HCTRecomb.lt.0.0) HCTRecomb = 0.0d0

      RETURN

      END


******************************************************************************

      DOUBLE PRECISION FUNCTION HeCTRecomb(nelem,ion,Te)

c
c     author: Chris Churchill
c
c     returns the helium charge exchange recombination rate coefficient
c     (cm^3/s) using the work of Kingdon + Ferland (1996,ApJS,106,205)

c     HeCTRecomb is returned in units cm^3/s

      implicit none
      integer           NTreg,ipIon,ion,nelem,ireg,iT
      double precision  Te,Tmin,Tmax,T4 
      double precision  a,b,c,d

      double precision  He1r
      common/CTHeRecomb/He1r(6,4,30,3)
      integer           NjkHe1r
      common/NkjHe1/    NjkHe1r(4,30)


      HeCTRecomb = 0.0d0

c     if off the element charge return; if this is the neutral ion, then
c     no CE, return

      IF (nelem.lt.2.OR.nelem.gt.30) RETURN
      IF (ion.eq.1) RETURN

c     note the translation to ion-1 !!!

      ipIon = ion - 1

c     for ions with j>4 we use the statistical results of Ferland,
c     Korista, & Verner (1997,ApJ,481,L115) Eq. 4; the constant 0.542E-9
c     is the reduced mass of the ion in terms of the electron mass for
c     hydrogen interactions; this is the limiting value for large mass
c     ions (i.e., iron); NOTE: they state that a factor of X uncertainty
c     in this constant yields a factor of uncertainty of X^n in the
c     resulting ionization fractions for Ipion=n+5 (n>=1)

      IF (ipIon.gt.4) then  ! use statistical charge transfer for ion > 4
       HeCTRecomb = 0.54d-9 * float(ipIon)
       RETURN
      END IF

      NTreg = NjkHe1r(ipIon,nelem)

c     if only 1 T region and the first fitting parameter is null, then
c     no data for this ion; hoping that this is because the rate is
c     negligible!

      IF ((NTreg.eq.1).AND.(He1r(1,ipIon,nelem,NTreg).eq.0.0)) RETURN

c     to stat within the allowed T range, this sets the T to be used in
c     the fitting function equal to the limiting Tmin or limiting Tmax;
c     this avoids setting the rate coefficient to zero while avoiding
c     extrapolation

c     find the index corresponding to the temperature range that
c     brackets te
      
      iT = 0
      IF (NTreg.gt.1) then
       DO 03 ireg=1,NTreg
        Tmin  = He1r(5,ipIon,nelem,ireg)
        Tmax  = He1r(6,ipIon,nelem,ireg)
        IF ((Te.ge.Tmin).AND.(Te.le.Tmax)) iT = ireg 
 03    CONTINUE
      ELSE
       iT = 1
      END IF

c     if we failed to assign iT, then bracket the extreme limits in T;
c     else grab the appropriate limits for this T

      IF (iT.eq.0) then
       Tmin  = He1r(5,ipIon,nelem,1)
       Tmax  = He1r(6,ipIon,nelem,NTreg)
      ELSE
       Tmin  = He1r(5,ipIon,nelem,iT)
       Tmax  = He1r(6,ipIon,nelem,iT)
      END IF

c     if T<Tmin then return; T>Tmax then grab the extreme upper limit
c     and conver to T4; else we are in a good T zone so convert T to T4

       IF (Te.lt.Tmin) RETURN
       IF (Te.gt.Tmax) then
        T4 = Tmax * 1.0d-4
        iT = NTreg
       ELSE
        T4 = Te * 1.0d-4
       END IF

c     compute the rate from the fitting function for this T

      a = He1r(1,ipIon,nelem,iT)
      b = He1r(2,ipIon,nelem,iT)
      c = He1r(3,ipIon,nelem,iT)
      d = He1r(4,ipIon,nelem,iT)

      HeCTRecomb = 1.0d-9 * a*(T4**b)*(1.0d0 + c*exp(d*T4))

      IF (HeCTRecomb.lt.0.0) HeCTRecomb = 0.0d0

      RETURN

      END


******************************************************************************

      DOUBLE PRECISION FUNCTION HCTIon(nelem,ion,Te)

c     returns the hydrogen charge exchange ionization rate
c     coefficient (cm^3/s) using the work of Kingdon + Ferland
c     (1996,ApJS,106,205)

*     ion is stage of ionization, 1 for atom
*     nelem is atomic number of element, 3 up to 30
*     Example:  O + H+ => O+ + H is HCTIon(8,1,T)

c     HCTRecomb is returned in units cm^3/s

      implicit none
      integer           ipIon,ion,nelem
      double precision  Te,Tmin,Tmax,T4 
      double precision  a,b,c,d,Bf

      double precision  H2i
      common/CTIon/     H2i(7,3,30)



c     set the null value

      HCTIon = 0.0d0

c     note that there is no translation from ion to ipIon

      ipIon = ion

c     if off the element chart or ionization chart  return

      IF (nelem.lt.3.OR.nelem.gt.30.OR.iPIon.gt.3) RETURN

c     if the first fitting parameter is null, then no data for this ion;
c     hoping that this is because the rate is negligible!

      IF (H2i(1,ipIon,nelem).eq.0.0) RETURN

c     to stat within the allowed T range, this sets the T to be used in
c     the fitting function equal to the limiting Tmin or limiting Tmax;
c     this avoids setting the rate coefficient to zero while avoiding
c     extrapolation

      Tmin  = H2i(5,ipIon,nelem)
      Tmax  = H2i(6,ipIon,nelem)

      IF (Te.lt.Tmin) RETURN
      IF (Te.gt.Tmax) then
       T4 = Tmax * 1.0d-4
      ELSE
       T4 = Te * 1.0d-4
      END IF

c     compute the rate from the fitting function; note that the
c     ionization is determmined from detailed balancing of the
c     recombination modulated by the Boltzmann factor

      a  = H2i(1,ipIon,nelem)
      b  = H2i(2,ipIon,nelem)
      c  = H2i(3,ipIon,nelem)
      d  = H2i(4,ipIon,nelem)
      Bf = H2i(7,ipIon,nelem)

      HCTIon = 1.0d-9 * a*(T4**b)*(1.0d0 + c*exp(d*T4)) * exp(-Bf/T4)

      IF (HCTIon.lt.0.0) HCTIon = 0.0d0

      RETURN

      END



******************************************************************************

      SUBROUTINE Hesetup

c
c     author: Chris Churchill
c
c     this routine fills the matrix NjkHe1r, which stores the number of
c     temperature regimes over which the He recombination rate
c     coefficients have fitting parameters

c     NjkHe1r is the logical ength of the 4th dimension in He1r(j,k)
c     matrix

      integer           NjkHe1r
      common/NkjHe1/    NjkHe1r(4,30)


c     first, initialize all ions to have a single temperature regime

      DO 11 j=1,4
       DO 09 k=1,30
        NjkHe1r(j,k) = 1
 09    CONTINUE
 11   CONTINUE

c     manually populate the special cases

      NjkHe1r(4,4)  = 2
      NjkHe1r(2,5)  = 2
      NjkHe1r(4,5)  = 2
      NjkHe1r(4,6)  = 3
      NjkHe1r(2,7)  = 2
      NjkHe1r(4,7)  = 2
      NjkHe1r(2,8)  = 2
      NjkHe1r(2,9)  = 3
      NjkHe1r(3,9)  = 2
      NjkHe1r(4,9)  = 2
      NjkHe1r(2,10) = 3
      NjkHe1r(3,10) = 3
      NjkHe1r(4,10) = 2
      NjkHe1r(2,11) = 2
      NjkHe1r(3,11) = 3
      NjkHe1r(3,12) = 2
      NjkHe1r(3,17) = 2
      NjkHe1r(2,19) = 3
      NjkHe1r(3,19) = 3
      NjkHe1r(4,19) = 2
      NjkHe1r(3,20) = 3
      NjkHe1r(4,22) = 2
      NjkHe1r(4,27) = 2
      NjkHe1r(3,28) = 2
      NjkHe1r(3,29) = 2

      RETURN
      END

