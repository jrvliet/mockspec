c     PACKAGE: Mockspec/anal
c     PROGRAM: sysanal
c     MODULE:  calc.f
c
c     DESCRIPTION
c     routines that compute the various absorption quantities; this is
c     the heart and soul of sysanal
c
c
c     this file contains:
c     SUBROUTINE finalcalc
c     SUBROUTINE ewcalc
c     SUBROUTINE drcalc
c     SUBROUTINE aodcalc
c     SUBROUTINE vmomcalc
c
c
c.........................................................................
c

      SUBROUTINE finalcalc

c     driver for the final analysis of the systems

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i

      include            'sysanal.com'


c     compute or recompute all quantities

      DO 03 i=1,norders
       CALL aodspec(i)
       CALL ewcalc(i)
       CALL drcalc(i)
       CALL aodcalc(i)
       CALL vmomcalc(i)
 03   CONTINUE

c     return

      RETURN
      END


c
c.........................................................................
c

      SUBROUTINE ewcalc(iorder)

c     Part 1. measure the equivalent widths and their uncertainties in
c     each master region defined by the 1st ion in the list, include in
c     the calculation only the subregions within each master region
c     (this treatment minimizes the errors by eliminating insignificant
c     flux decrements); also compute the flux weighted central
c     wavelength and lambda width; when appropriate in, compute the
c     equivalent width limit, which is designated by setting the
c     uncertainty to -1

c     Part2. compute the total equivalent width for the ion transition,
c     being sure to account for limits

c     all units in angstroms; observed values are stored, not rest-frame
c     values. (rest-frame values are output in the results files
c     generated from SYSANAL)

c
c.........................................................................
c

      implicit none
    
      include             'sysanal.h'

      integer             i,j,k,iorder
      double precision    ew0,ew1,ew2,sig,wmom2,dlam,flxnrm,
     &                    signrm,flxdec,ew0dn,ew0up,ctmp,csig

      include             'sysanal.com'

 
      Slev = N_sigma(iorder)

c     IORDER is passed in the call; this routine computes the values for
c     one ion transition per call

c     zero the total system EW and tota system uncerntainty in the EW
c     for this ion transition

       ewtot(iorder)    = 0.0d0
       ewsigtot(iorder) = 0.0d0

c     NLINES is the number of master regions, loop over the lines one by
c     one and compute the quantities

       DO 51 k=1,nlines

c     initialize working variables for the current order

        ew0  = 0.0d0
        ew1  = 0.0d0       
        ew2  = 0.0d0       
        sig  = 0.0d0
        dlam = 0.0d0 

c     Part 1. compute the equivalent width (zeroth moment of the flux
c     decrement); compute the centroid of the line (first moment of the
c     flux decrement); compute the uncertainty in the equivalent width;
c     the outer loop is the number of subregions in the master region,
c     the inner loop is over the pixels in the subregion; this treatment
c     avoids pixels with insignificant flux decrements 

        DO 52 i=1,nfind(iorder,k)
         DO 53 j=f_beg(k,iorder,i),f_end(k,iorder,i)
          flxnrm = flux(j,iorder)/cont(j,iorder)
          signrm = sigma(j,iorder)/cont(j,iorder)  
          IF (j.lt.ndata(iorder)) then
            dlam = wave(j+1,iorder)-wave(j,iorder)
          ELSE
            dlam = wave(j,iorder)-wave(j-1,iorder)
          END IF
          flxdec = 1.0d0 - flxnrm
          ew0    = ew0 + dlam * flxdec
          ew1    = ew1 + dlam * wave(j,iorder) * flxdec
          sig    = sig  + (dlam * signrm)**2
 53      CONTINUE
 52     CONTINUE

c     account for uncertainty in the continuum fit; this method employs
c     the approximation of Sembach & Savage (1992) in which the
c     uncertainty is 0.3 times the sigma spectrum (they adopt 0.4);
c     first we adjust the continuum downward and compute ew0dn and then
c     we adjust the continuum upward and compute ew0up; the "systematic"
c     uncertainty due to continuum placement is then given by csig =
c     0.5*abs(ew0up-ew0dn)

        ew0dn = 0.0d0  ! downward continuum EW
        DO 102 i=1,nfind(iorder,k)
         DO 103 j=f_beg(k,iorder,i),f_end(k,iorder,i)
          ctmp   = cont(j,iorder) - 0.3*abs(sigma(j,iorder))
          flxnrm = flux(j,iorder)/ctmp
          IF (j.lt.ndata(iorder)) then
            dlam = wave(j+1,iorder)-wave(j,iorder)
          ELSE
            dlam = wave(j,iorder)-wave(j-1,iorder)
          END IF
          flxdec = 1.0d0 - flxnrm
          ew0dn  = ew0dn + dlam * flxdec
 103     CONTINUE
 102    CONTINUE

        ew0up = 0.0d0  ! upward continuum EW
        DO 104 i=1,nfind(iorder,k)
         DO 105 j=f_beg(k,iorder,i),f_end(k,iorder,i)
          ctmp   = cont(j,iorder) + 0.3*abs(sigma(j,iorder))
          flxnrm = flux(j,iorder)/ctmp
          IF (j.lt.ndata(iorder)) then
            dlam = wave(j+1,iorder)-wave(j,iorder)
          ELSE
            dlam = wave(j,iorder)-wave(j-1,iorder)
          END IF
          flxdec = 1.0d0 - flxnrm
          ew0up  = ew0up + dlam * flxdec
 105     CONTINUE
 104    CONTINUE

c     take the half difference

        csig = 0.50d0*abs(ew0up-ew0dn)

c     store the values; the systematic and statistical uncertainies are
c     added in quadrature

        ew(k,iorder)    = ew0
        ewsig(k,iorder) = sqrt(abs(sig**2 + csig**2))

c     if we have data here then store the data; if the signifcance level
c     is less than SLEV (runtime input parameter), then replace with the
c     limit at the SLEV level

        IF (ew0.ne.0.0d0) then
         wbar(k,iorder ) = ew1/ew0
         siglevel(k,iorder) = ew(k,iorder)/ewsig(k,iorder)
         IF (siglevel(k,iorder).lt.Slev) then
          ew(k,iorder)    = Slev*ewsig(k,iorder)
          ewsig(k,iorder) = -1.0
          sf_flag(iorder,k) = .false. 
         END IF

c     if we do not have data, then zero WBAR and move on

        ELSE

         wbar(k,iorder ) = 0.0d0

        END IF


c     compute the width of the feature (the second moment of the flux
c     decrement), this requires WBAR; the looping is the same as above

       DO 54 i=1,nfind(iorder,k)
         DO 55 j=f_beg(k,iorder,i),f_end(k,iorder,i)  
          flxnrm = flux(j,iorder)/cont(j,iorder)
          IF (j.lt.ndata(iorder)) then
            dlam = wave(j+1,iorder)-wave(j,iorder)
          ELSE
            dlam = wave(j,iorder)-wave(j-1,iorder)
          END IF
          wmom2  = dlam *  (wave(j,iorder) - wbar(k,iorder))**2
          flxdec = 1.0d0 - flxnrm
          ew2    = ew2 + wmom2 * flxdec
 55      CONTINUE
 54     CONTINUE

c     store the value of the width

        width(k,iorder) = sqrt(ew2)


c     next master region on this transition (order)

 51    CONTINUE


c     Part 2. compute the total equivalent width from all features in
c     the ion transition; if the feature is insignficant (now flagged by
c     SIGEW=-1) do not increment EWTOT but do increment the uncertainty
c     EWSIGTOT based upon the fact that the EW is the SLEV limit


       ew0 = 0.0d0
       sig = 0.0d0

       DO 59 k=1,nlines
         IF (ewsig(k,iorder).ne.-1.0) then ! add them only if they are detected
           ew0 = ew0 + ew(k,iorder)
           sig = sig + ewsig(k,iorder)**2
         END IF
 59    CONTINUE

c     if we have a non-zero ew0, then store the totals if we have ew0=0,
c     then we have a limit, compute the total limit and store; in DO
c     loop 61, note that the SIG is incremented using the EW, which may
c     seem incorrect; however, the above logic was that when a limit is
c     found, the uncertainty is set to -1 and the EWlimits was stored in
c     EW.

        IF (ew0.ne.0.0d0) then
          ewtot(iorder)    = ew0
          ewsigtot(iorder) = sqrt(sig)
          sltot(iorder)    = ewtot(iorder)/ewsigtot(iorder)
        ELSE
          sig = 0.0d0  ! for good measure
          DO 61 k=1,nlines
            sig = sig + ew(k,iorder)**2  ! see above comment
 61       CONTINUE
          ewtot(iorder)    = sqrt(sig)
          ewsigtot(iorder) = -1.0
          sltot(iorder)    =  0.0d0
        END IF

c     return

      RETURN
      END


c
c.........................................................................
c

      SUBROUTINE drcalc(iorder)

c     Compute the ratios of the equivalent widths and their
c     uncertainties for each master region relative to the 1st ion in
c     the list 

c     Compute the ratio for the total equivalent width for the ion
c     transition

c     all units in angstroms; no redsdhift dependence

c
c.........................................................................
c
      implicit none
    
      include             'sysanal.h'

      integer             k,iorder
      double precision    diorder,dnorder

      include             'sysanal.com'


c     IORDER is passed in the call; this routine computes the values for
c     one ion transition per call

c     use the 1st ion transition in the list for the regions; compute EW
c     ratios and uncertainty; account for limits

      DO 51 k=1,nlines

       dr(k,iorder) = ew(k,1)/ew(k,iorder)

       IF (ewsig(k,iorder).ne.-1.0) THEN
        diorder      = ewsig(k,iorder)/ew(k,iorder)
        dnorder      = ewsig(k,1)/ew(k,1)
        drsig(k,iorder) = dr(k,iorder)*sqrt(diorder**2 + dnorder**2)
       ELSE
        drsig(k,iorder) = -1.0d0        
       END IF

 51   CONTINUE

c     compute the total system EW ratios; account for limits

       drtot(iorder) = ewtot(1)/ewtot(iorder)

       IF (ewsigtot(iorder).ne.-1.0) THEN
        diorder          = ewsigtot(iorder)/ewtot(iorder)
        dnorder          = ewsigtot(1)/ewtot(1)
        drsigtot(iorder) = drtot(iorder)*sqrt(diorder**2 + dnorder**2)
       ELSE
        drsigtot(iorder) = -1.0d0        
       END IF

c     return

      RETURN
      END


c
c.........................................................................
c

      SUBROUTINE aodcalc(iorder)

c     Compute the integrated optical depth and column density for ion
c     transition IORDER in the master regions defined by the first ion
c     transition in the list; however compute only over the subregions
c     in these master regions; take care with both lower and upper
c     limits
 
c     if there are no detections in any of the regions, then we need to
c     compute the limits; done in a separate loop

c
c.........................................................................
c


      implicit none
    
      include             'sysanal.h'

      integer             i,j,k,ik,pix1,npix,nsat,iorder
      logical             satflag,limflag
      double precision    dv,dlam,base10
      double precision    c10,dcd10,dcu10
      double precision    tsum,tsup,tsdn,csum,csup,csdn
      double precision    nresel,lambar,dlambar

      include             'sysanal.com'



      i = iorder

c     store the conversion from base e to base 10 for the errors

      base10 = 1.0d0/log(10.0d0)

c     set the number of master regions from the first ion transition
c     inthe list

c     zero the summation elements; set the limits flag low

      tsum = 0.0d0
      tsup = 0.0d0
      tsdn = 0.0d0
      csum = 0.0d0
      csup = 0.0d0
      csdn = 0.0d0

      limflag = .false.

c     loop over the master regions

      DO 31 k=1,nlines

c     set the saturation flags low; we will assume that two contiguous
c     saturate pixesl constitues an upper limit on the optical depth and
c     column density; this can happen for very narrow unresolved
c     satruration

      nsat       = 0
      satflag    = .false.      

c     if the current ion transition has no detection, then skip and use
c     errors

       IF (sf_flag(i,k)) THEN

c     loop over the subregions in the current master region and sum over
c     the pixels in the subregion; there is a bit of an inconvenience
c     here; the column density spectra and their uncertainties are in
c     log10 space; we must unfold that for the summations, thus the C10,
c     DCD10, and DCU10 variables

        DO 51 ik=1,nfind(i,k)
         DO 53 j=f_beg(k,i,ik),f_end(k,i,ik)     

c     check upper limits and saturation flag

          IF ((collim(j,i).eq.1).AND.(nsat.eq.0)) then
           nsat = j
          END IF
          IF ((collim(j,i).eq.1).AND.(nsat.eq.j-1)) then
           nsat    = 0
           satflag = .true.
          END IF

c     avoid summing insignificant pixels and bad pixels; in these pixels
c     of the AOD spectra, the tau and N are the limits and the
c     uncertainties are vanished

          IF ((collim(j,i).ne.-1).AND.(collim(j,i).ne.-2)) then
           dv    = 0.5*(vel(j+1,i)-vel(j-1,i))
           tsum  = tsum + dv*tau(j,i)
           tsup  = tsup + (dv*dutau(j,i))**2
           tsdn  = tsdn + (dv*ddtau(j,i))**2
           c10   = 10.0d0**col(j,i)
           dcd10 = c10*ddcol(j,i)/base10
           dcu10 = c10*ducol(j,i)/base10
           csum  = csum + dv*c10
           csup  = csup + (dv*dcu10)**2
           csdn  = csdn + (dv*dcd10)**2
          END IF

 53      CONTINUE
 51     CONTINUE

       END IF

c     next master region (LINEI)

 31   CONTINUE


c     at this point if tsum=0 and csum=0, then we need to compute the
c     limits

      IF (tsum.eq.0.0) then

c     compute the lower limits; in these pixels of the AOD spectra, the
c     tau and N are the limits and the uncertainties are vanished; thus
c     we treat the tau and col values as uncertainties to get the limits
c     also, we use nresel to scale the AOD columns per resel

       limflag = .true.
       npix    = 0
       nresel  = 0
       lambar  = 0.0d0
       dlambar = 0.0d0

       DO 61 k=1,nlines
        DO 71 ik=1,nfind(i,k)
         DO 73 j=f_beg(k,i,ik),f_end(k,i,ik)
          lambar  = lambar + wave(j,i)
          dlambar = dlambar + 0.5*(wave(j+1,i)-wave(j-1,i))
          dv      = 0.5*(vel(j+1,i)-vel(j-1,i))
          tsum    = tsum + (dv*tau(j,i))**2
          c10     = 10.0d0**col(j,i)
          csum    = csum + (dv*c10)**2
          npix    = npix + 1
 73      CONTINUE
 71     CONTINUE
 61    CONTINUE

       lambar  = lambar/float(npix)
       dlambar = dlambar/float(npix)
       nresel = (lambar/dlambar) * (2.35*profile)

      END IF

c     now we have to set the totals and track limits; we have two
c     logical flags; 
c        LIMFLAG=true means we have upper limits
c        SATFLAG=true means we have saturation and lower limits
c        both these low means we have a good measurement

c     all subregions for this master region are included in the sum; set
c     the totals for the master region

       IF ((.not.limflag).AND.(.not.satflag)) then

c     we have good data

        tautot(i)   = tsum
        dutautot(i) = sqrt(tsup)
        ddtautot(i) = sqrt(tsdn)
        coltot(i)   = log10(csum)
        ducoltot(i) = base10*sqrt(csup)/csum
        ddcoltot(i) = base10*sqrt(csdn)/csum

       END IF

       IF ((.not.limflag).AND.(satflag)) then

c     we have some saturated pixels; upper limits

        tautot(i)   = tsum
        dutautot(i) = 1.0d0
        ddtautot(i) = 0.0d0
        coltot(i)   = log10(csum)
        ducoltot(i) = 1.0d0
        ddcoltot(i) = 0.0d0

       END IF

c     we have no significant pixels; lower limits at Slev (Nsigma) we
c     just use the average for the tau limit and use nresel for the
c     columns

       IF (limflag) then

        tautot(i) = Slev*sqrt(tsum/float(npix))
        dutautot(i) = 0.0d0
        ddtautot(i) = -1.0d0
        coltot(i) = log10(Slev*sqrt(nresel*csum/float(npix)))
        ducoltot(i) = 0.0d0
        ddcoltot(i) = -1.0d0

       END IF

c     return

      RETURN
      END


c
c.........................................................................
c

      SUBROUTINE vmomcalc(iorder)

c     Compute the velocity moments from the flux decrements over the
c     master regions but using only pixels in the subregions for each
c     master region in order to reduce errors (much like the equivalent
c     width and the AOD comptations)

c     we make extensive use of the SL_FLAG boolean; all values are
c     already vanished in routine zeroall; so if SL_FLAG is low (no
c     features in the master region) then we skip the calculation of the
c     velocity moments; this returns zero for master regions with no
c     detections in the current ion transition

c
c.........................................................................
c

      implicit none
    
      include             'sysanal.h'

      integer             i,j,k,ik,iorder
      double precision    v0,v1,v2,v3,sig0,sig1,sig2,sig3
      double precision    dv1dI,dv2dI,dv3dI,v0tot,v1tot
      double precision    dv,flxdec,dvbar

      include             'sysanal.com'



c     IORDER is passed in the call; this routine computes the values for
c     one ion transition per call

c     use the 1st ion transition in the list for the regions

       i = iorder

       v0tot   = 0.0d0
       v1tot   = 0.0d0

c     Part 1. compute the first three velocity moments of each region;
c     use the regions of the first ion transition in the list

c     before embarking on the computations; be sure the current ion
c     transition has detection; if it does not, then we would end up
c     computing noise; return leaving the results vanished

       DO 51 k=1,nlines

c     initialize all values for the current order; it vital to do this
c     prior to deciding if we skip this master region because these
c     terms are used for summations

         dv    = 0.0d0 
         v0    = 0.0d0
         v1    = 0.0d0       
         v2    = 0.0d0       
         v3    = 0.0d0
         sig0  = 0.0d0
         sig1  = 0.0d0
         sig2  = 0.0d0
         sig3  = 0.0d0
         dv1dI = 0.0d0
         dv2dI = 0.0d0
         dv3dI = 0.0d0
        
c     only if there is a detection here; if not we keep the zeroed value

        IF (sf_flag(i,k)) then

c     compute zeroth and 1st velocity moments; the outer loop is the
c     number of subregions in the master region, the inner loop is over
c     the pixels in the subregion; this treatment avoids pixels with
c     insignificant flux decrements


         DO 52 ik=1,nfind(i,k)
          DO 53 j=f_beg(k,i,ik),f_end(k,i,ik)
           dv     = 0.5 * (vel(j+1,i)-vel(j-1,i))
           flxdec = 1.0d0 - flux(j,i)/cont(j,i)
           v0     = v0 + dv * flxdec
           v1     = v1 + dv * vel(j,i) * flxdec
 53       CONTINUE
 52      CONTINUE 

         vbar(k,i)    = v1/v0

c     compute the second and third moments (requires vbar); and the
c     uncertainty integrals; the looping is the same as above.

         DO 54 ik=1,nfind(i,k)
          DO 55 j=f_beg(k,i,ik),f_end(k,i,ik)  
           dv     = 0.5 * (vel(j+1,i)-vel(j-1,i))
           dvbar  = vel(j,i) - vbar(k,i)
           flxdec = 1.0d0 - flux(j,i)/cont(j,i)
           v2     = v2 + dv*flxdec*dvbar**2 
           v3     = v3 + dv*flxdec*dvbar**3
 55       CONTINUE
 54      CONTINUE

         vwidth(k,i) = sqrt(abs(v2/v0))
         vasym(k,i)  = sign(1.0d0,v3)*abs(v3/v0)**(1.0d0/3.0d0)  

c     compute the uncertainties, the derivation of which is found in
c     Sembach and Savage (1992, ApJS, 83, 147); however, their
c     expressions are incorrect, there are no summations in the partial
c     derivatives wrt to the flux, they are individual terms for each
c     pixel

         DO 56 ik=1,nfind(i,k)
          DO 57 j=f_beg(k,i,ik),f_end(k,i,ik)  
           dv    = 0.5 * (vel(j+1,i)-vel(j-1,i))
           dvbar = vel(j,i) - vbar(k,i)
           sig0  = - dv/cont(j,i)
           sig1  = - vel(j,i)*dv/cont(j,i)
           sig2  = - dv*dvbar**2 / cont(j,i)
           sig3  = - dv*dvbar**3 / cont(j,i)
           dv1dI = dv1dI + ((v0*sig1 - v1*sig0)*sigma(j,i)/v0**2)**2
           dv2dI = dv2dI + ((v0*sig2 - v1*sig0)*sigma(j,i)/v0**2)**2
           dv3dI = dv3dI + ((v0*sig3 - v1*sig0)*sigma(j,i)/v0**2)**2
 57       CONTINUE
 56      CONTINUE

         sigvbar(k,i)   = sqrt(dv1dI)
         sigvwidth(k,i) = sqrt(dv2dI)/(2.0d0*vwidth(k,i))
         sigvasym(k,i)  = sqrt(dv3dI)/(3.0d0*vasym(k,i)**2)

c     next feature on this transition (order); increment for the total
c     first moment; which is used in Part 2 to compute the total second
c     and third moments

         v0tot   = v0tot + v0
         v1tot   = v1tot + v1

        END IF

 51    CONTINUE

c     Part 2. compute the above properties over the full system; which
c     requires vtotbar; skip the whole mess of there are no detections
c     on this ion transition

      IF (v0tot.ne.0.0) then

       v0    = v0tot
       v1    = v1tot
       v2    = 0.0d0       
       v3    = 0.0d0
       dv1dI = 0.0d0
       dv2dI = 0.0d0
       dv3dI = 0.0d0

       vtotbar(i) = v1/v0

       DO 61 k=1,nlines

c     do the total sums over the detected subregions; the moments for
c     each ion transition are computed realtive to the vbar for that
c     transition

        IF (sf_flag(i,k)) then

         DO 62 ik=1,nfind(i,k)
          DO 63 j=f_beg(k,i,ik),f_end(k,i,ik)
           dv     = 0.5 * (vel(j+1,i)-vel(j-1,i))
           dvbar  = vel(j,i) - vtotbar(i)
           flxdec = 1.0d0 - flux(j,i)/cont(j,i)
           v2     = v2 + dv*flxdec*dvbar**2 
           v3     = v3 + dv*flxdec*dvbar**3
 63       CONTINUE
 62      CONTINUE

        END IF

c     next feature on this transition (order)

 61    CONTINUE

       vtotwidth(i)    = sqrt(abs(v2/v0))
       vtotasym(i)     = sign(1.0d0,v3)*abs(v3/v0)**(1.0d0/3.0d0) 

c     compute the uncertainties over all master regions

       DO 64 k=1,nlines

c     do the total sums over the detected subregions

        IF (sf_flag(i,k)) then

         DO 65 ik=1,nfind(i,k)
          DO 66 j=f_beg(k,i,ik),f_end(k,i,ik)  
           dv     = 0.5 * (vel(j+1,i)-vel(j-1,i))
           dvbar  = vel(j,i) - vtotbar(i)
           sig0   = - dv/cont(j,i)
           sig1   = - vel(j,i)*dv/cont(j,i)
           sig2   = - dv*dvbar**2 / cont(j,i)
           sig3   = - dv*dvbar**3 / cont(j,i)
           dv1dI = dv1dI + ((v0*sig1 - v1*sig0)*sigma(j,i)/v0**2)**2
           dv2dI = dv2dI + ((v0*sig2 - v1*sig0)*sigma(j,i)/v0**2)**2
           dv3dI = dv3dI + ((v0*sig3 - v1*sig0)*sigma(j,i)/v0**2)**2
 66       CONTINUE
 65      CONTINUE

         END IF
 64    CONTINUE

       sigvtotbar(i)   = sqrt(dv1dI)
       sigvtotwidth(i) = sqrt(dv2dI)/(2.0d0*vtotwidth(i))
       sigvtotasym(i)  = sqrt(dv3dI)/(3.0d0*vtotasym(i)**2)

      END IF
     
c     return

      RETURN
      END

