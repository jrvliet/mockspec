c
c.........................................................................
c

      SUBROUTINE photoxsecs

c     populates the arrays of the photoionization cross sections for
c     each shell s of each ion k,j; we also populates each their energy
c     arrays, which must start with the threshold energy and then
c     populate the 2nd derivatives later required for interpolation for
c     the photoionization rate integrals (see function PH_rate)


c     the partial cross sections are obtained via a call to a modified
c     version of the code, phfit2.f, originally written by D. Verner; it
c     has now been thoroughly commented and modified; mods include
c     changing it from a subroutine call to a function call (thus the
c     argument returned, the cross section, has been removed from the
c     call parameter list); the module is now called phfit2-mod.f, see
c     comments therein for details (especially about the bug I found!)

c     the shell numbers are indexed (by s) as follows
c     s=1 1s (l=0)
c     s=2 2s (l=0)
c     s=3 2p (l=1)
c     s=4 3s (l=0)
c     s=5 3p (l=1)
c     s=6 3d (l=2)
c     s=7 4s (l=0)

c     we assume k such that the species numbers k go as 1=H, 2=He, etc,
c     and the ionization stages go as 1=neutral, 2=single, 3=double
c     ionized etc.; the ion is assumed to be in its ground state; our
c     formalism is to index each ion by (k,j), whereas Verner's notation
c     for his phfit2 routine indexes each element as (nz,ne); the
c     translation is
c
c     OUR'S                  VERNER'S
c     k                      nz       atomic number
c     nek                    ne       number of electrons bound to ion
c     s                      is       shell number
c
c     such that Verner's indices in terms of our indices are
c     nz = k               = atomic number
c     ne = nek = k - j + 1 = number of bound electrons
c     is = s

c     note this routine is caled only if doPH is set high and is called
c     after the SED arrays are populated; if we are using the UVB or the
c     UVB+stars SEDs, then the energy grid for the cross sections will
c     be based upon the UVB energy grid; if stars only are being used,
c     then the energy grid will be based upon stars only

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           i,j,k,kk,nek,s
      include           'rates.com'
      include           'com-photo.com'
      include           'com-phxsecs.com'

c     we need the number of shells in each isoelectronic sequence

c     we need the total number of populated shells, so include the ntot

      integer           ntot
      COMMON/ntot/      ntot(30)


c     loop over included species k and their ionization stages j;
c     populate the cross section arrays by energy, k, j, and shell s

      DO kk=1,Nspecies
       k = kidx(kk)
       DO j=1,k
        nek = k - j + 1
        DO s=1,ntot(nek)
         CALL popxkjs(k,j,s)
         CALL popd2xkjsdE2(k,j,s)
        ENDDO
       ENDDO  ! next j
      ENDDO   ! next kk
 
      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE popxkjs(k,j,s)

c     populate the photoioniation cross section arrays 

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           i,ie,iie,j,k,nek,s
      double precision  E
      double precision  phfit2,thresholdE
      include           'rates.com'
      include           'com-modes.com'
      include           'com-photo.com'
      include           'com-phxsecs.com'


c     number of electrons for ground state of ion k,j

      nek = k - j + 1

c     store the first energy element as the threshold energy

      E                   = thresholdE(k,nek,s)
      logEsigkjs(1,k,j,s) = log10(E)
      logsigkjs(1,k,j,s)  = log10(phfit2(k,nek,s,E)) - 18.0d0

      iie = 1

      IF (uvbflag) then  ! also if both ubvflag and sb99flag

c     populate to the highest energies

       DO ie=1,NEuvb
        IF (logEuvb(ie).gt.logEsigkjs(1,k,j,s)) then
         iie                   = iie + 1
         E                     = 10.0d0**logEuvb(ie) 
         logEsigkjs(iie,k,j,s) = logEuvb(ie) 
         logsigkjs(iie,k,j,s)  = log10(phfit2(k,nek,s,E)) - 18.0d0
        END IF
       ENDDO      
       NEsigkjs(k,j,s) = iie   ! number of data points

      ELSE  ! if sb99flag only

       DO ie=1,NEsb99
        IF (logEsb99(ie).gt.logEsigkjs(1,k,j,s)) then
         iie = iie + 1
         E                     = 10.0d0**logEsb99(ie)
         logEsigkjs(iie,k,j,s) = logEsb99(ie)
         logsigkjs(iie,k,j,s)  = log10(phfit2(k,nek,s,E)) - 18.0d0
        END IF
       ENDDO      
       NEsigkjs(k,j,s) = iie   ! number of data points

      END IF


      RETURN
      END


c
c.........................................................................
c

      SUBROUTINE popd2xkjsdE2(k,j,s)

c     modified Numerical Recipe routine, computes 1D 2nd derivatives of
c     array Y with respect to array X

c     this routine called once; so we translate the input arrays into
c     the canned routine (which is slower, but only done at the
c     beginning of the program so does not affect run time during
c     iterations on cells and subscells)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           k,j,s,i,n
      double precision  YP1,YPN,SIG,P,QN,UN
      double precision  X(NEmax),Y(NEmax),Y2(NEmax),U(NEmax)
      include           'com-phxsecs.com'


c     translate to the generalized notation of the canned routine

      N   = NEsigkjs(k,j,s) 
      YP1 = 1.0d30
      YPN = 1.0d30

      DO i=1,N
       x(i) = logEsigkjs(i,k,j,s)
       y(i) = logsigkjs(i,k,j,s)
      ENDDO

c     begin canned routine

      IF (YP1.GT..99D30) THEN
        Y2(1) = 0.0d0
        U(1)  = 0.0d0
      ELSE
        Y2(1) = -0.50d0
        U(1)  = (3.0d0/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF

      DO 11 I=2,N-1
        SIG = (X(I)-X(I-1))/(X(I+1)-X(I-1))
        IF (SIG.le.0.0) GOTO 99
        P     = SIG*Y2(I-1)+2.0d0
        Y2(I) = (SIG-1.0d0)/P
        U(I)  = (6.0d0*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *         /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE

      IF (YPN.GT..99D30) THEN
        QN = 0.0d0
        UN = 0.0d0
      ELSE
        QN = 0.50d0
        UN = (3.0d0/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF

      Y2(N) = (UN-QN*U(N-1))/(QN*Y2(N-1)+1.0d0)

      DO 12 i=N-1,1,-1
        Y2(i) = Y2(i)*Y2(i+1)+U(i)
12    CONTINUE

c     populate the derivatives from the canned routine

      DO i=1,N
       d2sigkjsdE2(i,k,j,s) = Y2(i)
      ENDDO

      RETURN

 99   WRITE(6,*) 'ERROR(xsecderivs): X data ill-posed'
      WRITE(6,*) 'index of trouble data (k,j,s,i) =',k,j,s,I
      WRITE(6,*) 'X data must be in ascending order'
      STOP

      END

     
c
c.........................................................................
c

      DOUBLE PRECISION FUNCTION phxkjs(E)

c     performs 1D interpolation of of the photoionization cross section
c     in log10 space

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           k,j,s,i
      integer           ilo,ihi
      double precision  a,b,hh,E,logE,y,y0,y1
      include           'com-phxsecs.com'
      include           'com-photo.com'

c     tranlate globals to locals

      k = kp
      j = jp
      s = sp

      logE = log10(E)

c     if this is the edge, just grab the cross section 

      IF (logE.eq.logEsigkjs(1,k,j,s)) then
       phxkjs = logsigkjs(1,k,j,s)
       RETURN
      ENDIF

c     otherwise, search for the index

      ilo = 1
      ihi = NEsigkjs(k,j,s) 


1     IF (ihi-ilo.GT.1) THEN
        i=(ihi+ilo)/2
        IF (logEsigkjs(i,k,j,s).GT.logE) THEN
          ihi=i
        ELSE
          ilo=i
        ENDIF
      GOTO 1
      ENDIF

      hh = logEsigkjs(ihi,k,j,s) - logEsigkjs(ilo,k,j,s)

      IF (hh.EQ.0.) THEN
       WRITE(6,*) 'ERROR(phxsec): logE cannot be bracketed due to'
       WRITE(6,*) 'bad logEsigkjs input at indices ilo,ihi= ',ilo,ihi
       STOP
      END IF

      a      = (logEsigkjs(ihi,k,j,s)-logE)/HH
      b      = (logE-logEsigkjs(ilo,k,j,s))/HH

      phxkjs = a*logsigkjs(ilo,k,j,s)+b*logsigkjs(ihi,k,j,s)+
     *         ((a**3-a)*d2sigkjsdE2(ilo,k,j,s)+
     *         (b**3-b)*d2sigkjsdE2(ihi,k,j,s))*(hh**2)/6.0d0

c     suppress ringing near discontinuities by enforcing consistency
c     with the overall shape; if ringing is occuring, use linear
c     interpolation

      y0 = min(logsigkjs(ilo,k,j,s),logsigkjs(ihi,k,j,s))
      y1 = max(logsigkjs(ilo,k,j,s),logsigkjs(ihi,k,j,s))

      IF ((phxkjs.lt.y0).OR.(phxkjs.gt.y1)) phxkjs = y0 + (y1-y0)*b

      RETURN
      END



c
