c     PACKAGE: Mockspec/anal
c     PROGRAM: sysanal
c     MODULE:  spectra.f
c
c     DESCRIPTION
c     routines that create the equivalent width spectra, EW(lambda), and
c     the apparent optical depth spectra, N(v); the uncertainty spectra
c     are also created
c
c
c     this file contains:
c     SUBROUTINE ewspec
c     SUBROUTINE aodspec
c
c
c.........................................................................
c

      SUBROUTINE ewspec(iorder)

c
c     We construct the equivalent width and sigma spectra using
c     Schneider etal 1993 (ApJS 87, 45); this spectrum will be searched
c     for finding absorption features.

c     CAVEAT: if your sigma spectrum is under (over) predicting your
c     uncertainties in your flux values then this routine will over
c     (under) find significant features; be sure your sigma spectrum is
c     sound

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i,j,k,kk,iorder
      double precision   isf(mxlin),sigisf,flxnrm,signrm,signrmk2
      double precision   bige,bigs,dw,norm,prob

      include           'sysanal.com'



c     initializing

      DO 09 i=1,ndata(iorder)
       ewpix(i,iorder)    = 0.0d0
       ewsigpix(i,iorder) = 0.0d0
 09   CONTINUE

c     loop over the innner regions of the spectrum that do not require
c     special handling of the ISF; inside this loop, when we are on the
c     first or last pixel, we treat the order ends

      DO 11 i=Jnot+1,ndata(iorder)-(Jnot+1)

c     construct the ISF and normalize; the ISF is constructed around
c     each pixel because it changes as a function of wavelength

       sigisf       = profile * wave(i,iorder)
       isf(Jnot+1)  = 1.0d0
       norm         = 1.0d0
       DO 17 j=1,Jnot
        dw   = (wave(i-j,iorder) - wave(i,iorder))/sigisf
        isf(Jnot+1-j) = exp(-0.5d0*(dw**2))
        dw   = (wave(i+j,iorder) - wave(i,iorder))/sigisf
        isf(Jnot+1+j) = exp(-0.5d0*(dw**2))
        norm = norm + isf(Jnot+1+j) + isf(Jnot+1-j)
 17    CONTINUE
       isf(Jnot+1) = isf(Jnot+1)/norm
       DO 18 j=1,Jnot
        isf(Jnot+1-j) = isf(Jnot+1-j)/norm
        isf(Jnot+1+j) = isf(Jnot+1+j)/norm
 18    CONTINUE
        
c     the "convolution" index k here is the equivalent to their
c     i+j-Jnot, which is corrected here for an error by replacing 
c     j -> j-1

       bige = 0.0d0
       bigs = 0.0d0
       prob = 0.0d0
       dw   = 0.5 * abs(wave(i+1,iorder)-wave(i-1,iorder))

       DO 13 j=1,M
        k      = i + (j-1) - Jnot
        flxnrm = flux(k,iorder)/cont(k,iorder)
        signrm = sigma(k,iorder)/cont(k,iorder)
        bige   = bige + isf(j) * (flxnrm-1.0d0)  
        bigs   = bigs + (isf(j)*signrm)**2
        prob   = prob + isf(j)**2
 13    CONTINUE

       ewpix(i,iorder)    = dw * bige / prob
       ewsigpix(i,iorder) = dw * sqrt(bigs) / prob 

c     if we are on the first pixel for which the isf is complete...  do
c     the front edge now; what we are assuming here is that the
c     continuum stretches beyond the order edge with a value of unity;
c     we "borrow" the first dw, isf, and prob values

       IF (i.eq.Jnot+1) then

         DO 21 kk=1,Jnot

          bige = 0.0d0
          bigs = 0.0d0

          DO 23 j=1,M
           k        = kk + (j-1) - Jnot
           flxnrm   = flux(k,iorder)/cont(k,iorder)
           signrm   = sigma(k,iorder)/cont(k,iorder)
           signrmk2 = sigma(kk,iorder)/cont(kk,iorder)
           IF (k.lt.1) then
             bigs = bigs + (isf(j)*signrmk2)**2
           ELSE
             bige = bige + isf(j) * (flxnrm-1.0d0)  
             bigs = bigs + (isf(j)*signrm)**2
           END IF 
 23       CONTINUE

          ewpix(kk,iorder)    = dw * bige / prob
          ewsigpix(kk,iorder) = dw * sqrt(bigs) / prob 

 21      CONTINUE

       END IF

c     if we are on the last pixel for which the isf is complete...  do
c     the back edge now; what we are assuming here is that the continuum
c     stretches beyond the order edge with a value of unity; we "borrow"
c     the last dw, isf, and prob values

       IF (i.eq.ndata(iorder)-(Jnot+1)) then

         DO 31 kk=ndata(iorder)-Jnot,ndata(iorder)

          bige = 0.0d0
          bigs = 0.0d0

          DO 33 j=1,M
           k        = kk + (j-1) - Jnot
           flxnrm   = flux(k,iorder)/cont(k,iorder)
           signrm   = sigma(k,iorder)/cont(k,iorder)
           signrmk2 = sigma(kk,iorder)/cont(kk,iorder)
           IF (k.gt.ndata(iorder)) then
             bigs = bigs + (isf(j)*signrmk2)**2
           ELSE
             bige = bige + isf(j) * (flxnrm-1.0d0)  
             bigs = bigs + (isf(j)*signrm)**2
           END IF 
 33      CONTINUE

          ewpix(kk,iorder)    = dw * bige / prob
          ewsigpix(kk,iorder) = dw * sqrt(bigs) / prob 

 31     CONTINUE

       END IF


 11   CONTINUE

c     return

      RETURN
      END


c.........................................................................
c

      SUBROUTINE aodspec(iorder)

c
c     Compute the optical depth spectrum for ion transition IORDER in
c     the master regions defined by the first ion transition in the list
c
c.........................................................................
c


      implicit none
    
      include             'const.dek'
      include             'sysanal.h'

      integer             j,i,iorder,pix1,npix
      double precision    flxnrm,signrm,con
      double precision    zb,ze,z1,z2,base10

      include             'sysanal.com'



c     store the conversion from base e to base 10 for the errors

      base10 = 1.0d0/log(10.0d0)

c     compute the optical depth constant, requires physical constants
c     input from const.dek; the factor of 1.E13 converts the speed of
c     light to km/s (1.E5) and the wavelength to centimeters (1.E8)

      con = 1.0d13*me*c/(pi*e*e*fval(iorder)*lambda0(iorder))

c     compute tau in each pixel uncertainty; account for limits; we have
c     vanished tau of we are outside a detection region; we use the
c     uncertainty if we are saturated

c     loop over master regions from first ion transition in the list

      DO 31 j=1,nlines

c     we compute the spectra over the full master region projected onto
c     the current ion transition; grab the pixel limits

       pix1 = 0
       npix = 0

       zb = wave(f_beg(j,1,1),1)/lambda0(1)
       ze = wave(f_end(j,1,1),1)/lambda0(1)
       DO 09 i=1,ndata(iorder)-1
        z1 = wave(i,iorder)/lambda0(iorder)
        z2 = wave(i+1,iorder)/lambda0(iorder)
        IF ((pix1.eq.0).AND.(zb.lt.z1)) pix1 = i   ! low-z edge effect (Q&D)
        IF ((npix.eq.0).AND.(ze.gt.z2)) npix = i   ! high-z edge effect (Q&D)
        IF ((zb.ge.z1).AND.(zb.le.z2)) pix1 = i
        IF ((ze.ge.z1).AND.(ze.le.z2)) npix = i
 09    CONTINUE

       IF (pix1.eq.0) then
        WRITE(6,*) 'ERROR(aodspec): pix1=0; cannot locate region'
        STOP
       ENDIF
       IF (npix.eq.0) then
        WRITE(6,*) 'ERROR(aodspec): npix=0; cannot locate region'
        STOP
       ENDIF

c     loop over the pixels for this region

       DO 51 i=pix1,npix

c     bad pixels are assumed to have sigma=-1.0; this will biff the
c     optical depth calculation; if pad pix, retain the value

        flxnrm = flux(i,iorder)/cont(i,iorder)
        IF (sigma(i,iorder).ne.-1.0d0) then
         signrm = sigma(i,iorder)/cont(i,iorder)  
        ELSE
         signrm = sigma(i,iorder)
        END IF

c     trap bad pixels, assumed to have SIGNRM=-1.0

        IF (signrm.eq.-1.0d0) then
         tau(i,iorder)    = 0.0d0
         col(i,iorder)    = 0.0d0
         collim(i,iorder) = 0
         GOTO 41
        END IF

c     trap saturation? upper limits are set; skip other traps

        IF ((flxnrm.lt.signrm).OR.(flxnrm.le.0.0)) then
         tau(i,iorder)    = -log(signrm)
         col(i,iorder)    = con*tau(i,iorder)
         collim(i,iorder) = 1 
         GOTO 41
        END IF

c     trap insignificant (consistent with continuum); lower limits are
c     set; skip other traps

        IF (flxnrm+signrm.gt.1.0d0) then
         tau(i,iorder)   = -log(1.0d0-signrm)
         col(i,iorder)   = con*tau(i,iorder)
         collim(i,iorder) = -1 
         GOTO 41
        END IF

c     "light's green, trap is clean"; if we are here: pixel value is
c     solid

        tau(i,iorder)   = -log(flxnrm)
        dutau(i,iorder) = -0.5*log(flxnrm-signrm)
        ddtau(i,iorder) = -0.5*log(flxnrm+signrm)
        col(i,iorder)   = con*tau(i,iorder)
        ducol(i,iorder) = con*dutau(i,iorder)
        ddcol(i,iorder) = con*ddtau(i,iorder)

c     set the uncertainties and convert the column density to log10;
c     trap pad pixels

 41     IF (signrm.ne.-1.0d0) then
          ducol(i,iorder) = base10*ducol(i,iorder)/col(i,iorder)
          ddcol(i,iorder) = base10*ddcol(i,iorder)/col(i,iorder)
          col(i,iorder)   = log10(con*tau(i,iorder))
        ELSE
          ducol(i,iorder) = 0.0d0
          ddcol(i,iorder) = 0.0d0
          col(i,iorder)   = 0.0d0
        END IF

 51    CONTINUE


 31   CONTINUE

c     return

      RETURN
      END

c     eof
