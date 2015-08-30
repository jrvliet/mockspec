c.........................................................................
c

      SUBROUTINE timescales(T,Zmet,kdo,jdo)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           j,k,m,Ndo
      integer           kdo(Imax),jdo(Imax)
      double precision  yr2sec,kbolt
      parameter         (yr2sec = 3.15569d7,
     &                   kbolt  = 1.380658d-16)
      double precision  coolfunc,tauphoto,taurecomb,taucoll,acoll
      double precision  T,Zmet,tau,Qg,Qgdot,Zg
      include           'rates.com'


c     gas energy (Qg), and energy loss rate (Qgdot); the Qgdot
c     computation assumes collisional 
c     **** not sure about Zg=Zmet at this stage ****

      Zg     = Zmet
      Qg     = 1.50d0*(eden+nions)*kbolt*T
      Qgdot  = nions*nions*coolfunc(T,Zg)

c     cooling time scale (assuming CIE, which is not correct!)

      tau      = Qg/Qgdot
      tau_cool = tau/yr2sec
      tau_cool = log10(tau_cool)

c     loop over output ions

      DO m=1,noutfiles

c     obtain the ion indices

       k = kdo(m)
       j = jdo(m)

c     if k=0, this is the electron density, skip function calls

       IF (k.eq.0) then
        tau_ph(m)   = 1.0d-50
        tau_rec(m)  = 1.0d-50
        tau_coll(m) = 1.0d-50
        GOTO 01 
       ENDIF

c     function calls to compute the timescales

       tau_ph(m)   = tauphoto(k,j)/yr2sec    ! photo
       tau_rec(m)  = taurecomb(k,j)/yr2sec   ! recombination
       tau_coll(m) = taucoll(k,j)/yr2sec     ! collisional

c     convert to log10 for output/communication

 01    tau_ph(m)   = log10(tau_ph(m))
       tau_rec(m)  = log10(tau_rec(m))
       tau_coll(m) = log10(tau_coll(m))

      ENDDO

C      WRITE(6,*) nions,eden,Zg,T,coolfunc(T,Zg),Qg,Qgdot,tau_cool


      RETURN
      END


c.........................................................................
c

      DOUBLE PRECISION FUNCTION coolfunc(T,Zg)

c     compute the specific cooling function for CIE; Maller & Bullock
c     (2004) cooling function approximation (with corrections)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      double precision Lam23,T,cfunc,Zg
      double precision Tb,Tm,Tr,alpha
      double precision third,twothird

      parameter       (third    = 1.0d0/3.0d0,
     &                 twothird = 2.0d0/3.0d0,
     &                 Tm       = 1.50000d5,
     &                 Tr       = 1.50000d4)

c     Tr = temperature boundary for radiative cooling
c     Tm = temperature transition to metal cooling
c     Tb = temperature boundary for bremsstrahlung cooling
c     alpha = slope in metal cooling regime

      Tb    = 1.0d6 + 1.5d7*(Zg**twothird)
      alpha = 1.0d0 + third*log(Zg)

      IF ((T.ge.Tr).AND.(T.le.Tm)) then
       coolfunc = 1.0d-23*12.0d0*((T/Tr)**alpha)
       RETURN
      ENDIF

      IF ((T.gt.Tm).AND.(T.le.Tb)) then
       Lam23     = 12.0d0*((Tm/Tr)**alpha)
       coolfunc  = 1.0d-23*Lam23*(Tm/T)
       RETURN
      ENDIF

      IF (T.gt.Tb) then
       Lam23    = 12.0d0*((Tm/Tr)**alpha) * (Tm/Tb)
       coolfunc = 1.0d-23*Lam23*((T/Tb)**third)
       RETURN
      ENDIF

c     temporary fix for T<Tr

      IF (T.lt.Tr) then
       coolfunc = 1.0d-26
       RETURN
      ENDIF

c     death trap
      WRITE(6,*) 'ERROR(coolfunc): temperature out of range'
      WRITE(6,*) 'T = ',T 
      STOP

      END


c.........................................................................
c

      DOUBLE PRECISION FUNCTION tauphoto(k,j)

c     returns the photoionization time scale

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           j,k
      double precision  efold,nrat
      include           'rates.com'


      nrat  = nkj(k,j-1)/nkj(k,j)

c     check if this is the neutral ion

      IF (j.eq.1) then
       efold = R_ph(k,j)
       IF (efold.ne.0.0) then
        tauphoto = 1.0d0/efold
       ELSE
        tauphoto = 3.15569d-43   ! "null"
       ENDIF
       RETURN
      ENDIF
      
c     metal ions (most likely)

      efold = abs(R_ph(k,j) - nrat*R_ph(k,j-1))
      IF (efold.ne.0.0) then
       tauphoto = 1.0d0/efold
      ELSE
       tauphoto = 3.15569d-43   ! "null"
      ENDIF

      RETURN
      END

c.........................................................................
c

      DOUBLE PRECISION FUNCTION taurecomb(k,j)

c     returns the recombination time scale 

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           j,k
      double precision  efold,nrat
      include           'rates.com'


      nrat  = nkj(k,j+1)/nkj(k,j)

c     check if this is the neutral ion

      IF (j.eq.1) then
       efold = nrat*eden*alpha_rec(k,j)
       IF (efold.ne.0.0) then
        taurecomb = 1.0d0/efold
       ELSE
        taurecomb = 3.15569d-43   ! "null"
       ENDIF
       RETURN
      ENDIF

c     metal ions (most likely)

      efold = abs(eden*(alpha_rec(k,j-1)-nrat*alpha_rec(k,j)))
      IF (efold.ne.0.0) then
       taurecomb = 1.0d0/efold
      ELSE
        taurecomb = 3.15569d-43   ! "null"
      ENDIF

      RETURN
      END


c.........................................................................
c

      DOUBLE PRECISION FUNCTION taucoll(k,j)

c     returns the collisional ionization time scale

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           j,k
      double precision  efold,nrmj,nrpj,acollj,acollmj
      include           'rates.com'


      nrmj    = nkj(k,j-1)/nkj(k,j)
      nrpj    = nkj(k,j+1)/nkj(k,j)
      acollj  = alpha_cdi(k,j) + alpha_cea(k,j)
      acollmj = alpha_cdi(k,j-1) + alpha_cea(k,j-1)

c     check if this is the neutral ion

      IF (j.eq.1) then
       efold = abs(eden*acollj)
C       efold = abs(eden*(acollj - alpha_rec(k,j)))
       IF (efold.ne.0.0) then
        taucoll = 1.0d0/efold
       ELSE
        taucoll = 3.15569d-43   ! "null"
       ENDIF
       RETURN
      ENDIF

c     metal ions (most likely)

       efold = abs(eden*(acollj - nrmj*acollmj))
C       efold = abs(eden*(acollj + alpha_rec(k,j-1) -
C     &               nrmj*acollmj - nrpj*alpha_rec(k,j)))
      IF (efold.ne.0.0) then
       taucoll = 1.0d0/efold
      ELSE
       taucoll = 3.15569d-43   ! "null"
      ENDIF

      RETURN
      END
