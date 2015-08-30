c
c
c     IMPORTANT NOTES ON ENERGY UNITS OF SEDS and CROSS SECTIONS: the
c     SED is used for calculating the photoionization rate coefficients;
c     these calculations involve integrating the product of the
c     photoionization cross section and the mean intensity of the SED
c     over photon energies greater than than the threshold for
c     ionization; this integral (performed in funcs-photo.f) is done
c     with the mean intensity in units [cm^-1 s^-1 Str^-1], effectively
c     the mean intensity density per unit energy interval.  Since energy
c     units do not appear in the mean intensity density per unit energy
c     interval, the integration can be perfomed in any energy unit;
c     i.e., photon energies in eV, or in erg, etc.  We use eV,
c     throughout; the cross sections are in units [cm^2].
c
c     Thus, for all arrays involving spectral energy distributions (UVB,
c     Stars, etc), we first convert the SED to the mean intensity
c     density per unit energy interval and store the energy for each
c     point in eV.

c
c.........................................................................
c

      SUBROUTINE mkUVBofz(z)

c     this routine called only if user has set the boolean "uvbflag" to
c     a hight state

c     interpolate the UVB spectra at various redshifts to obtain the
c     spectrum at the desired redshift, z

c     we do this by interpolating in the z direction at each energy
c     point

c     we return log(JEuvb), not log(JNUuvb), so that the integrals can
c     be done over energy and not frequencies

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           ie,iz,n
      double precision  z,yp1,ypN,x,y
      double precision  xx(UVBmaxpix),yy(UVBmaxpix),y2(UVBmaxpix)
      include           'com-photo.com'

c     set up the natural cubic spline

      yp1 = 1.0d30
      ypN = yp1

c     we interpolate at each energy point in the grid; there NEuvb
c     energy pixels and n=UVBzpix redshift pixels

      DO 09 ie=1,NEuvb

       DO 11 iz=1,UVBzpix
        xx(iz) = zuvb(iz)
        yy(iz) = logJNUuvb(iz,ie)
        y2(iz) = 0.0d0
 11    CONTINUE

       n = UVBzpix  ! convert to a local variable for spline calls
       CALL spline(xx,yy,n,yp1,ypN,y2)
       CALL splint(xx,yy,y2,n,z,y)

c     convert logJnu via Jnu*dnu=J_E*dE, dnu/dE=1/h to the the mean
c     intensity density per unit energy interval [cm^-1 s^-1 Str^-1]

       logJEuvb(ie) = y - log10(h) 

 09   CONTINUE


      RETURN
      END


c
c.........................................................................
c

      SUBROUTINE mkSB99spop(spop_age,spop_Zmet,spop_Msol,rstar)

c     this routine called only if user has set the boolean "sb99flag" to
c     a hight state

c     this routine takes the grid of stellar population spectra and
c     interpolates to obtain the spectrum for the desired age,
c     metallicity, and mass 

c     the input spectra are in units erg s^-2 Ang^-1 (i.e., they are
c     luminosity density)

c     currently, there is no mass grid because the liminosity scales
c     with the number of stars and the number of stars scales with mass,
c     so we can do a simple scaling law; the grid currently uses a
c     stellar population with 10^3 solar masses

c     we first interpolate across ages for each metallicity in the grid
c     to obtain the spectra at the desired age at each grid metallicity

c     we then interpolate across metallicity at the desired age

c     finally, we convert the spectrum to the mean intensity per unit
c     energy, JE, erg cm^-2 str^-1, as we did for the UVB; this is done
c     in two steps, first we convert to the flux per unit Ang (this
c     requires the inverse square law), and then we convert to JE


c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           ipix,iZ,iage,n
      double precision  spop_age,spop_Zmet,spop_Msol,rstar
      double precision  yp1,ypN,y
      double precision  logLageZ,Rcm,logFlam
      dimension         logLageZ(SB99NZ,SB99maxpix)
      double precision  logLuse(SB99maxpix),swap(SB99maxpix)
      double precision  xx(SB99maxAZ),yy(SB99maxAZ),y2(SB99maxAZ)
      include           'com-photo.com'



c     set up the natural cubic spline

      yp1 = 1.0d30
      ypN = yp1

c     first, interpolate each metallicity to obtain the spectrum at the
c     desired age; we interpolate at each wavelength point in the grid
c     across age; there are SB99NZ metallicity pixels, NEsb99 energy
c     pixels, and SB99age age pixels

      DO 01 iZ=1,SB99NZ
       DO 03 ipix=1,NEsb99

        DO 05 iage=1,SB99Nage
         xx(iage) = Agesb99(iage)
         yy(iage) = logLlamsb99(iage,iZ,ipix)
         y2(iage) = 0.0d0
 05     CONTINUE

        n = SB99Nage  ! convert to a local variable for spline calls
        CALL spline(xx,yy,n,yp1,ypN,y2)
        CALL splint(xx,yy,y2,n,spop_age,y)

c     store the luminosity density for the desired age 

        logLageZ(iZ,ipix) = y 

 03    CONTINUE
 01   CONTINUE

c     second interpolate the luminosity density at the desire age to
c     obtain the luminosity density for both the desired age and the
c     desired metallicity

c     we perform a linear interpolation for now until we get more
c     metallicities on the grid

      DO 13 ipix=1,NEsb99
       logLuse(ipix) = logLageZ(1,ipix) 
     &         + ((spop_Zmet-Zsolsb99(1))*logLageZ(2,ipix)   
     &         -  (spop_Zmet-Zsolsb99(1))*logLageZ(1,ipix))    
     &         / (Zsolsb99(2)-Zsolsb99(1))
 13   CONTINUE


c     third, convert the luminosity density per unit wavelength to the
c     flux density per unit wavelength, and convert the flux density per
c     unit wavelength to the mean intensity per unit energy

      Rcm = rstar * 1.0d3 * pc2cm

      DO 15 ipix=1,NEsb99
       logFlam = logLuse(ipix) - log10(4.0d0*pi) - 2.0d0*log10(Rcm)
       logJEsb99(ipix) = logFlam + log10(h*clight) + 8.0d0
     &      - 2.0d0*(logEsb99(ipix)+log10(erg2eV)) - log10(4.0d0*pi)
 15   CONTINUE

c     last, the SB99 spectra are stored in ascending wavelenght order,
c     whic hmeans that they are in descending energy order; we now
c     invert the order

c     swap the E array

      n = 0
      DO 21 ipix=NEsb99,1,-1
       n = n + 1
       swap(n) = logEsb99(ipix)
 21   CONTINUE
      n = 0
      DO 22 ipix=1,NEsb99
       logEsb99(ipix) = swap(ipix)
 22   CONTINUE

c     swap the JE array

      n = 0
      DO 23 ipix=NEsb99,1,-1
       n = n + 1
       swap(n) = logJEsb99(ipix)
 23   CONTINUE
      n = 0
      DO 24 ipix=1,NEsb99
       logJEsb99(ipix) = swap(ipix)
 24   CONTINUE


      RETURN
      END


c     storing the computation of log10(energy) to log10(lambda)
c     --------------------------------------------------------
c        lambda  = log10(h*clight) - logEsb99(ipix) + 8.0d0
c     &                 - log10(erg2eV)
c        y   = logFlam - log10(4.0d0*pi)
c        WRITE(77,*) logEsb99(ipix),logJEsb99(ipix),y
c     --------------------------------------------------------

c
c.........................................................................
c

      SUBROUTINE mkSED

c     add the UVB and SB99 ionizing spectra 

c     this routine called only if user has set the boolean "doPH" to a
c     hight state; also, if doPH=.false. then the routine PH_rate is not
c     called
 
c     ASSUMPTION: the UVB spectrum covers more energy range than the SB99
c     spectrum

c     we have multiple conditions to consider

c     (1) we are using a combined UVB plus stellar contribution
c     (2) we are using only the UVB contribution
c     (3) we are using only the stellar contribution

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           ie,n
      double precision  Elosb99,Ehisb99,Elouvb,Ehiuvb,logE
      double precision  y,yp1,ypN,yuvb,ysb99
      double precision  xx(SEDmaxpix),yy(SEDmaxpix),y2(SEDmaxpix)
      include           'com-photo.com'
      include           'com-modes.com'

c     set up the natural cubic spline

      yp1 = 1.0d30
      ypN = yp1

c     ..............................................................
c     BEGIN OPTION 1: combining both the UVB and the stellar spectra
c     ..............................................................

      IF ((uvbflag).AND.(sb99flag)) then  ! begin OPTION 1

c     first copy the UVB spectrum to the SED spectrum for all energies

       NEsed = NEuvb  ! assumption: UVB covers broader energy range
       DO 15 ie=1,NEsed
        logEsed(ie)  = logEuvb(ie)
        logJEsed(ie) = logJEuvb(ie)
 15    CONTINUE

c     now, we interpolate the SB99 spectrum to the energy grid of the
c     UVB spectrum in regions of energy overlap and add them (replacing
c     the copied values above in the overlap regions)

       DO 19 ie=1,NEsb99
        xx(ie) = logEsb99(ie)
        yy(ie) = logJEsb99(ie)
        y2(ie) = 0.0d0
 19    CONTINUE
       n = NEsb99
       CALL spline(xx,yy,n,yp1,ypN,y2)
       Elosb99 = logEsb99(1)
       Ehisb99 = logEsb99(NEsb99)
       DO 21 ie=1,NEsed
        logE = logEsed(ie)
        IF (logE.ge.Elosb99.AND.logE.lt.Ehisb99) then
         CALL splint(xx,yy,y2,n,logE,y)
         yuvb         = 10.0d0**logJEuvb(ie)
         ysb99        = 10.0d0**y
         y            = yuvb + ysb99
         logJEsed(ie) = log10(y)
        END IF
 21    CONTINUE

      ENDIF  

c     END OPTION 1

c     ...............................................
c     BEGIN OPTION 2: using only the UVB contribution
c     ...............................................

      IF ((uvbflag).AND.(.not.sb99flag)) then  ! begin OPTION 2

c     simply copy the UVB spectrum to the SED spectrum for all energies

       NEsed = NEuvb  
       DO 25 ie=1,NEsed
        logEsed(ie)  = logEuvb(ie)
        logJEsed(ie) = logJEuvb(ie)
 25    CONTINUE

      ENDIF

c     END OPTION 2

c     ...............................................
c     BEGIN OPTION 3: using only the UVB contribution
c     ...............................................

      IF ((.not.uvbflag).AND.(sb99flag)) then  ! begin OPTION 3

c     simply copy the SB99 spectrum to the SED spectrum for all energies

       NEsed = NEsb99
       DO 27 ie=1,NEsed
        logEsed(ie)  = logEsb99(ie)
        logJEsed(ie) = logJEsb99(ie)
 27    CONTINUE

      ENDIF

c     END OPTION 3

c     .......................................
c     finalize setup of the ionizing spectrum
c     .......................................

c     obtain the global 2nd derivatives in the energy direction of
c     logJEuvb for future interpolating along the energy direction

      CALL SEDderivs

c     this sets the energy buffer so that we have constraining data for
c     interpolating the cross section 10-11 data points past the energy
c     cut off of the Jnu spectrum (makes the interpolation more robust)

      logEmax = logEsed(NEsed)
      Ebuff   = 10.0d0*(logEsed(NEsed) - logEuvb(NEsed-1))

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE SEDderivs

c     computes 1D 2nd derivatives of array Y with respect to array X

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           i,n
      double precision  U(SEDmaxpix),UN,QN,SIG,P
      include           'com-photo.com'
      include           'com-modes.com'

C       xx(ie) = logEsed(ie)
C       yy(ie) = logJEsed(ie)
C       d2JEsed(ie) = y2(ie)

      N = NEsed

      d2JEsed(1) = 0.0d0
      U(1)       = 0.0d0
      QN         = 0.0d0
      UN         = 0.0d0

      DO 11 i=2,N-1
        SIG = (logEsed(i)-logEsed(i-1))/(logEsed(i+1)-logEsed(i-1))
        IF (SIG.le.0.0) GOTO 99
        P          = SIG*d2JEsed(i-1)+2.0d0
        d2JEsed(i) = (SIG-1.0d0)/P
        U(I)       = (6.0d0*((logJEsed(I+1)-logJEsed(I))/(logEsed(I+1)-
     *               logEsed(I))-(logJEsed(I)-logJEsed(I-1))
     *               /(logEsed(I)-logEsed(I-1))) / 
     *               (logEsed(I+1)-logEsed(I-1))-SIG*U(I-1))/P
11    CONTINUE

      DO 12 i=N-1,1,-1
        d2JEsed(i)=d2JEsed(i)*d2JEsed(i+1) + U(i)
12    CONTINUE

      d2JEsed(N)=(UN-QN*U(N-1))/(QN*d2JEsed(N-1)+1.0d0)

      RETURN

 99   WRITE(6,*) 'ERROR(spline): logEsed data ill-posed index I=',I
      WRITE(6,*) 'X data must be in ascending order'
      STOP

      END

