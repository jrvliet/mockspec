c
c.........................................................................
c

      DOUBLE PRECISION FUNCTION PH_Rate(k,j)

c     this routine has two purposes:
c
c     (1) compute the partial photoionization rate for species k in
c     ionization stage j (by partial we mean for each individual shell);
c     store these in the array R_phs(s), where "s" is the shell
c     number (see below); the partial rates will be used to compute the
c     Auger rates

c     (2) compute the total photoionization rate (summed over all shells)
c     for single electron ejection (no Auger electrons); return this
c     rate as PH_Rate

c     the partial photoionization rate coefficient is determined by
c     integrating the ionizing spectrum modulated by the partial
c     photoionization cross section (units cm^2) over energy from the
c     threshold energy for the shell to the maximum energy of the
c     ionizing spectrum
c
c     the shell numbers are indexed (by s) as follows
c     s=1 1s (l=0)
c     s=2 2s (l=0)
c     s=3 2p (l=1)
c     s=4 3s (l=0)
c     s=5 3p (l=1)
c     s=6 3d (l=2)
c     s=7 4s (l=0)
c
c     the partial cross sections are populated in routine photoxsecs

c     the integral for obtaining the rates is 
c
c      integral = 4*pi * integral{sigkj(s,E)*J(E)*dE/E}   E=E0->Emax
c

c     the integration is performed in routine rate_integral; we perform
c     the integration in parts to speed up the convergence

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           i,j,k,nek,s,nout,Nsplit
      double precision  logElo,logEhi,dlogE,Elo,E,pow
      double precision  Eth(NShmax)
      double precision  alphakj,aphsum
      double precision  phfit2,thresholdE
      include           'rates.com'
      include           'com-photo.com'
      include           'com-phxsecs.com'
      include           'com-auger.com'

c     we need the total number of populated shells, so include the ntot

      integer           ntot
      COMMON/ntot/      ntot(30)

c     ...........
c     SET UP SHOP
c     ...........

c     determine the number of bound electrons, nek

      kp   = k
      jp   = j
      nek  = k - j + 1

c     .....................
c     BEGIN THE INTEGRATION
c     .....................

c     integrate over cross section and ionizing spectral distribution to
c     obtain the photoionization rate coefficient one shell at a atime

      DO 01 s=1,ntot(nek)

       sp     = s
       aphsum = 0.0d0

c     we must trap the special cases for the KI, CaI and CaII electron
c     configurations (4s1, 4s2) since the shell s=6 (3d) is empty, so
c     null and advance to the next shell

       IF (s.eq.6) then
        IF ((k.eq.19.AND.nek.eq.19).OR.    ! KI    (4s1)
     &      (k.eq.20.AND.nek.eq.19).OR.    ! CaII  (4s1)
     &      (k.eq.20.AND.nek.eq.20)) then  ! CaI   (4s2)
         R_phs(s) = 0.0d0
         GOTO 01  ! advance to next shell
        ELSE
         GOTO 02  ! carry on
        END IF
       END IF

c     set up the energy grid for the partial cross section calculation;
c     Emin is the current shell energy and is the lower energy limit
c     over which the integration for the rate is performed; Emax is
c     global and has been returned from routine getJEofz (the maximum
c     energy of the ionizing spectrum)

 02    logEmin = logEsigkjs(1,k,j,s)

c     if the shell energy is greater than the maximum energy of the SED
c     array, then then the contribution from this shell is null, skip to
c     next shell (the shell energies decline as shell number incerases)

       IF (logEmin.ge.logEsed(NEsed)) GOTO 01

c      if the shell energy is less than SED min energy then we have a
c      situation where ionization from this shell can still occur ans so
c      it will contribute to the photoionization rate, but we are
c      missing the lowest energy photons from the SED array, we fix this
c      by resetting the minimum energy of the integration to the lowest
c      energy in the SED array; note however, that this represent some
c      missed contribution to the photoionization rate in principle, but
c      not from a computational standpoint (i.e., in the real world the
c      SED energy does not cut off to absolute zero.. however, as can be
c      seen with stellar only SEDS, the energy drops so rapidly that
c      this treatment is still an excellent approximation)

       IF (logEmin.lt.logEsed(1)) logEmin = logEsed(1)

c     good to go, set up the energy grid

       dlogE  = 0.01
       logElo = logEmin
       logEhi = logEmax + Ebuff
       NExsec = (logEhi-logElo)/dlogE + 1

c     integrate beyond the current shell edge one energy decade at a
c     time; this splitting of the integral speeds up convergence because
c     it is difficult to numerically integrate a function over several
c     (4-6!) orders of magnitude in energy; variable Nsplit determines
c     the number of full energy decades to loop over

       logEmax = logEsed(NEsed)
       dlogE   = logEmax - logEmin
       Nsplit  = int(dlogE)
       IF (Nsplit.eq.0) Nsplit = 1
       DO 04 i=1,Nsplit
        logEmax = logEmin + 1.0d0
        CALL rate_integral(alphakj)
        aphsum = aphsum + alphakj
        logEmin = logEmax
 04    CONTINUE

c     now integrate over the final energy decade (this energy interval
c     will not be a full decade- and the contribution is usually
c     negligible, though that depends upon the ioninzing SED)

       logEmin = logEmax
       logEmax = logEsed(NEsed)
       CALL rate_integral(alphakj)
       aphsum = aphsum + alphakj

c     store the rate for this shell

       R_phs(s) = aphsum

 01   CONTINUE  ! increment to next shell


c     ..........................................................
c     COMPUTE PHOTOIONIZATION RATE FOR SINGLE ELECTRON EJECTIONS
c     ..........................................................

c     the photoionization rate for single electron ejection is the sum
c     over all the shells with each shell weighted by the probability of
c     yielding a single electron (no Auger electrons)

      PH_Rate = 0.0d0
      DO 05 s=1,ntot(nek)
       PH_Rate = PH_Rate + W_Auger(k,j,s,1)*R_phs(s)
 05   CONTINUE

c     return

      RETURN

      END


c
c.........................................................................
c

      DOUBLE PRECISION FUNCTION thresholdE(k,nek,s)

c     find the threshold E; this is due to some freaking esoteric logic
c     in Verner's phfti2 code (which came uncommented!).  the need for
c     it is explained in the introduction to his 1996 paper; however, it
c     is not clear how it is implemented in his code

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           k,nek,s,is
      double precision  phfit2
      double precision  E,x
      include           'com-photo.com'

c     we need the threshold energies for each shell, so we include the
c     ph1 common block from Verner's PH1 block data; we also need the
c     inner shells so we include Verner's NINN block data
      
      double precision  ph1
      COMMON/ph1/       ph1(6,30,30,7) 
      integer           ninn
      common/ninn/      ninn(30)       ! inner shell



c     obtain the threshold for this shell as employed by Verner

      is = s
      E  = ph1(1,k,nek,is) 
      x  = phfit2(k,nek,is,E)

      IF (x.gt.0.0) then
       thresholdE = E
      ELSE
       is = ninn(nek)
       thresholdE = ph1(1,k,nek,is) 
      END IF

      RETURN
      END



c
c.........................................................................
c

      SUBROUTINE rate_integral(alpha)

c     integrate!  we integrate the function FUNK from E1 to E2

c     the integrand FUNK is the product of the photoionization cross
c     section and the mean intensity of the ionizing spectrum (SED)

c     we use the Numerical Recipe integrator called QROMB
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      external          funk
      double precision  funk
      double precision  E1,E2,alpha
      include           'com-photo.com'

      E1 = 10.0d0**logEmin

      E2 = 10.0d0**logEmax

      CALL qromb(funk,E1,E2,alpha)

      alpha = 4.0d0*pi*alpha

      RETURN
      END


c
c.........................................................................
c

      DOUBLE PRECISION FUNCTION funk(E)

c     the integrand, which is the product of the photoionization cross
c     section and the mean intensity of the ionizing spectrum (SED)

c     since we have stored the cross sections and the SED in log10 units
c     we exponentiate; recall that the cross section and the SED are
c     stored in arrays with dependent variable log10(E) at intervals of
c     E; thus, we require interpolation to obtain their values at
c     arbitrary E

c     the function FlogJEsed(E) performs cubic spline interpolation in
c     order to return log10 of the ionizing SED at arbitrary energy E

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      double precision  E,sigE,JE
      double precision  phxkjs,FlogJEsed
      include           'com-photo.com'

      
      sigE = 10.0d0**phxkjs(E)
      JE   = 10.0d0**FlogJEsed(E)

      funk = sigE*JE/E

      RETURN
      END


c
c.........................................................................
c

      DOUBLE PRECISION FUNCTION FlogJEsed(E)

c     this function returns the value of log(JEsed) at arbitrary energy
c     E (eV) for the integrating routine by spline interpolation

c     this is a converted Numerical Recipe (splint.f)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      IMPLICIT DOUBLE PRECISION (A-H,O)
      include           'rates.h'
      include           'com-photo.com'

      X = log10(E)

      KLO=1
      KHI=NEsed

1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF (logEsed(K).GT.X) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

      DH = logEsed(KHI) - logEsed(KLO)

      IF (DH.EQ.0.) THEN
       WRITE(6,*) 'ERROR(FlogJE): E cannot be bracketed due to'
       WRITE(6,*) 'bad logEsed input at indices KLO,KHI= ',KLO,KHI
       WRITE(6,*) 'E(KLO),JE(KLO) = ',logEsed(KLO),logJEsed(KLO)
       WRITE(6,*) 'E(KHI),JE(KHI) = ',logEsed(KHI),logJEsed(KHI)
       STOP
      END IF

      A = (logEsed(KHI)-X)/DH
      B = (X-logEsed(KLO))/DH
      Y = A*logJEsed(KLO)+B*logJEsed(KHI) +
     *   ((A**3-A)*d2JEsed(KLO)+(B**3-B)*d2JEsed(KHI))*(DH**2)/6.0d0

c     suppress ringing at ionization edges by enforcing consistency with
c     the SED shape; if ringing is occuring, use linear interpolation

      y0 = min(logJEsed(KLO),logJEsed(KHI))
      y1 = max(logJEsed(KLO),logJEsed(KHI))

      IF ((Y.lt.y0).OR.(Y.gt.y1)) Y = y0 + (y1-y0)*B

      FlogJEsed = Y

      RETURN
      END


