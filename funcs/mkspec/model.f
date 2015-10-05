c     PACKAGE: Mockspec
c     PROGRAM: specsynth
c     MODULE:  model.f
c
c
c     DESCRIPTION
c     these routines setup the information needed to make the model
c     spectrum and make the smoothed model spectrum (the spectrum
c     returned is not yet convolved with the ISF)
c
c
c     this file includes:
c     SUBROUTINE setatomic
c     SUBROUTINE initspectrum
c     SUBROUTINE doabslines
c     SUBROUTINE dolymanlimit
c     SUBROUTINE droplines
c
c
c..............................................................................
c

      SUBROUTINE setatomic(j)

c     here we assige the atomic constants to each line; I know this
c     seems archaic since the lines are all for one transition, but the
c     left over logic from the original general specsynth was that each
c     line could be a different transition; so as to not to need to
c     rewrite all the model.f routines; this was the compromise

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include           'specsynth.h'
      include           'const.dek'
      integer           i,j
      character*80      trani


c     assign an atomic constant to each line

      DO 13 i=1,lines
       lambda0(i) = lam0(j)
       con1(i)    = 1.0d+5*f0(j)*lam0(j)**2 * sqrt(pi)*e*e/(me*c*c)
       con2(i)    = 1.0e-8*gamma0(j)*lam0(j)**2 / (4.0d0*pi*c)
  13  CONTINUE

      write(*,*) con1(1), con2(1)
      RETURN

      END


c
c..............................................................................
c  
      SUBROUTINE initspectrum

c     as it says, initializes the arrays that house the spectrum

c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include           'specsynth.h'
      integer           pixi

c     given the wavelength range and the pixel sampling rate, determine
c     the number of pixels, "ndata", in the spectrum

      ndata = int((wave_max-wave_min)/dwave) + 1

c     sanity check; are the number of pixels greater than the allocated
c     physical sizes of the arrays; the parameter "maxpix" is set in the
c     "specsynth.h" file and can be modified therein

      IF (ndata.gt.maxpix) then
        WRITE(SCREEN,*) 'ERROR(initspectrum): NDATA > MAXPIX'
        STOP
      END IF

c     fill the wavelength pixel values and initialize the spectrum as
c     unity continuum in array "wrkflux"

      DO 11 pixi=1,ndata
       lambda(pixi) = wave_min + (pixi-1)*dwave
       wrkflx(pixi) = 1.0d0
 11   CONTINUE
  
c     return

      RETURN

      END


c
c..............................................................................
c  
      SUBROUTINE doabslines
  
c     used Voigt profile formalism and background source radiative
c     transfer (no source function) to compute the absorption lines

c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h'
      include           'const.dek'
      integer           linei,pixi
      double precision  w0,b1,b2,b3,w,x,y,u,v,tcon

c     some helpful definitions
c     b1  = column density
c     b2  = rest frame central wavelength
c     b3  = rest frame Doppler width
c     y   = natural broadening term in doppler units
c     x   = wavelength difference in doppler units
c     w0  = rest frame wavelength for transition
c     w   = rest frame wavelength for observed wavelength
c     tau = optical depth

c     main loop: loop over the individual lines one at a time and
c     generate the Voigt profile contribution in each pixel (inner loop)

      DO 05 linei=1,lines
  
        IF (INDEX(ion_name(linei),'LLS').ne.0) GOTO 05

c     initialize Voigt profile parameters; the parameters must be rest
c     frame quantities even though the line may be redshifted

        w0   = lambda0(linei)
        b1   = nline(linei)
        b2   = w0 * (1.0d0+zline(linei))/(1.0d0+zabs)
        b3   = w0 * abs(bline(linei)) /ckms
        y    = con2(linei) / b3
        tcon = con1(linei) * (b1/b3)

        write(*,*) w0, b1, b2, b3, y, tcon
c     inner loop: loop over the pixels and perform the radiative
c     transfer

        DO 11 pixi=1,ndata
          w = lambda(pixi)/(1.0+zabs)
          x = (w-b2)/b3 
          CALL voigt(x,y,u,v)
          wrkflx(pixi) = wrkflx(pixi) * exp(-tcon*u)
 11     CONTINUE

 05   CONTINUE

c     return

      RETURN

      END

  
c
c..............................................................................
c  

      SUBROUTINE dolymanlimit
  
c     if we are modelling the limit limit, then we are also modeling the
c     Lyman series (done in routine "doabslines"); this routine performs
c     the radiative transfer for the ionization edge of HI using the
c     semi-classical Kramer's ionization cross section 

c     I WANT TO REPLACE KRAMER WITH THE FULL CROSS SECTION  May 2011

c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h'
      include           'const.dek'
      integer           linei,pixi
      double precision  w0,w,tau,NHI,zHI


c     main loop: loop over the line one at a time; check if the line is
c     the break (designated "LLS"); we work only on the "LLS" "line"

      DO 03 linei=1,lines

        IF (INDEX(ion_name(linei),'LLS').eq.0) GOTO 03   ! if "T" skip

c       define the necessary quanitites for the cross section

        NHI = 1.0e13 * nline(linei)  
        zHI = zline(linei) 
        w0  = lambda0(linei)   ! this is the LL according to the atomic data
        tau = 0.0d0

        DO 11 pixi=1,ndata
          w  = lambda(pixi)/(1.0+zHI)
          IF (w.le.w0) then    ! then we are have ionizing photons! 
            tau = NHI * 6.30d-18 * (w/w0)**3
            wrkflx(pixi) = wrkflx(pixi) * exp(-tau)
          END IF
          IF ((w.gt.w0).AND.(w.lt.912.645)) then ! quick fix smooth the break
            tau = NHI * 6.3d-18 
            wrkflx(pixi) = wrkflx(pixi) * exp(-tau)
          END IF
 11     CONTINUE

 03   CONTINUE

c     return

      RETURN

      END
  

c
c..............................................................................
c
      SUBROUTINE droplines

c     drop lines that are not in the wavelength range; we do this to
c     make the spectrum generation more efficient in the case that there
c     are many mnay lines but that they are not in the desired
c     wavelength region specified by the user

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'specsynth.h'

      integer           i,k,dropped
      double precision  w


c     communicate the number of line originall yin the list

C      WRITE(6,'(a,i5)') ' Lines in input list         ',lines

c     initialize counters

      k       = 1
      dropped = 0

c     loop over lines and remove those outside the wavelength range

      DO 21  while (k.le.lines)

        w = lambda0(k)*(1.0d0+zline(k))
        IF ((w.lt.wave_min).OR.(w.gt.wave_max)) then
         lines = lines - 1
         dropped = dropped + 1
         DO 11 i=k,lines
           ion_name(i) = ion_name(i+1)
           lambda0(i)  = lambda0(i+1)
           con1(i)     = con1(i+1)
           con2(i)     = con2(i+1)
           zline(i)    = zline(i+1)
           nline(i)    = nline(i+1)
           bline(i)    = bline(i+1)  
 11      CONTINUE
        ELSE
         k = k + 1      ! count this as a included line
        END IF

 21   CONTINUE

c     communicate the number of lines dropped and the number of lines to
c     be included; I find this helps in the case that the use made bad
c     specifications for the desired spectrum or something was not
c     initialized properly with the instrument set up

C      WRITE(6,'(a,i5)') ' Lines out of spectral range ',dropped
C      WRITE(6,'(a,i5)') ' Lines being processed       ',lines

c     return

      RETURN

      END 

