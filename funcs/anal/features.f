c     PACKAGE: Mockspec/anal
c     PROGRAM: sysanal
c     MODULE:  features.f
c
c     DESCRIPTION
c     routines that objectively determine the regions in the spectra
c     where significant absorption features are located
c
c
c     this file contains:
c     SUBROUTINE initregs
c     SUBROUTINE spurious    ! not active
c     SUBROUTINE delregion   ! not active
c     INTEGER FUNCTION sysfeat
c     INTEGER FUNCTION subfeat
c
c
c.........................................................................
c

      SUBROUTINE initregs

c     Initial computations on the spectra prior to editing

c     This routine is called once at the start of the program; 

c     (1) creates the EW spectra

c     (2) automatigically finds the master regions, 

c     (3) automagically finds the subregions for the higher order ion
c         transitions, 

C     NOT ACTIVE
c     (4) determines which of the master regions are spurious (have no
c         features in higher order ion transitions), 

c     (5) computes the apparent optical depth (AOD) spectra, and 

c     (6) determines the absorber redshift 

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i,linei,sysfeat,subfeat,action
      logical            flag

      include            'sysanal.com'



c     compute the equivalent width spectra

      DO 01 i=1,norders
       CALL ewspec(i)
 01   CONTINUE

c     find the master regions; fill the R_BEG and R_END arrays that
c     store the starting and ending pixels of the system regions

      nlines = sysfeat()
      IF (nlines.eq.0) then
C       WRITE(6,*) ' No features detected'
       RETURN
      END IF

C     DEBUG code
C      DO linei=1,nlines
C      WRITE(6,*) linei,wave(f_beg(linei,1,1),1),wave(f_end(linei,1,1),1)
C      ENDDO

c     find the regions in higher order ion transitions, if there are
c     any?

      IF (norders.gt.1) then 
       DO 02 i=2,norders
        DO 03 linei=1,nlines
         nfind(i,linei) = subfeat(i,linei)
 03     CONTINUE
 02    CONTINUE
      END IF

C     this is commented out because it does not apply for the Mock
C     spctra, which are clean

c     checking for spurious master regions can be performed only if there
c     are higher order ion transitions

C      IF (norders.gt.1) then
C       flag   = .false.
C       action = 1
C       DO 04 linei=1,nlines
C        CALL spurious(linei,action,flag)
C 04    CONTINUE
C       IF (nlines.eq.0) then
C        WRITE(6,*) ' Features were spurious'
C        RETURN
C       END IF
C      END IF

c     deteremine the system redshift; this must be called after the
c     aodspec are loaded

      DO 05 i=1,norders
       CALL aodspec(i)
 05   CONTINUE

      CALL getzabs

c     return

      RETURN
      END


c
c.........................................................................
c

      INTEGER FUNCTION sysfeat()

c     This routine called once from routine INITREGS

c     This function returns the number of significant absorption
c     features, which we call the master regions, and the starting and
c     ending pixels, stored in f_beg and f_end, of each feature.  The
c     method uses the equivalent width spectrum, not the simple flux
c     decrement.  Though it can be used for emission features to, this
c     subroutine is designed to scope out absorption features only

c     A feature us found when one pixel has an EW in its pixel that is
c     Nsig beyond the EWSIG uncertainty in that pixel, i.e., EW/EWSIG <
c     -Nsig.  Then the spectrum is scanned backwards to find the blue
c     edge of the feature, which is defined when the ratio EW/EWSIG = -
c     1.  This is the location at about which the flux decrement becomes
c     consistent with the continuum.  Then the spectrum is set back to
c     the initial pixel of the detection and scanned forward for
c     EW/EWSIG = -1.  The pixel locations of the lower and upper bounds
c     are stored in F_BEG and F_END.  Then the pixel number is set to
c     the upper pixel + 1 of the just found feature and the scanning for
c     another Nsig detection is searched for.

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i,ipix,lpix,upix,npix
      integer            iorder,linei
      double precision   Nsig,ratio

      include            'sysanal.com'



      sysfeat = 0
      linei    = 0

c     set the ion transistion to the first one in the list; IORDER=1

      iorder   = 1

c     set the significance level

      Nsig     = N_sigma(iorder)

c     scan the full data chunk

      ipix     = 1
      npix     = ndata(iorder)

c     we scan the spectrum; ipix is the pixel meeting the detection
c     criterion; then we scan blueward and then redward looking for the
c     extremities of the feature; we use a while loop because we need to
c     leap frog the feature once we have it defined, thus we do not
c     increment ipix uniformily

      DO i=1,npix
       ratio = ewpix(i,iorder)/ewsigpix(i,iorder)
      ENDDO

      DO 11 WHILE (ipix.lt.npix)

       ratio = ewpix(ipix,iorder)/ewsigpix(ipix,iorder)

c     enter scan feature block if we have a significant equivalent
c     width; increment the feature index and then scan for the feature
c     extremities

       IF (ratio.le.-Nsig) then 

         linei = linei + 1
         lpix  = ipix
         upix  = ipix

c     search blueward for the beginning of the feature

         DO 13 WHILE ((ratio.lt.-1.0).AND.(lpix.gt.1))
           lpix  = lpix - 1 
           ratio = ewpix(lpix,iorder)/ewsigpix(lpix,iorder)  
 13      CONTINUE

c     reset the ratio and pixel to the scanning start point

         ratio = ewpix(ipix,iorder)/ewsigpix(ipix,iorder) 

c     search redward for the ending of the feature

         DO 14 WHILE ((ratio.lt.-1.0).AND.(upix.lt.npix))
           upix  = upix + 1 
           ratio = ewpix(upix,iorder)/ewsigpix(upix,iorder)  
 14      CONTINUE

c     front edge effects; warning

         IF (lpix.eq.1) then 
           WRITE(6,*) ' WARNING: feature',linei,' artificially ',
     &                'terminated start of ',order(iorder)(1:20)
           WRITE(1,*) ' WARNING: feature',linei,' artificially ',
     &                'terminated start of ',order(iorder)(1:20)
         END IF

c     back edge effects; warning

         IF (upix.eq.npix) then 
           WRITE(6,*) ' WARNING: feature',linei,' artificially ',
     &                'terminated end of ',order(iorder)(1:20)
           WRITE(1,*) ' WARNING: feature',linei,' artificially ',
     &                'terminated end of ',order(iorder)(1:20)
         END IF

c     got one- store the pixel indices and set the SL_FLAG high

         nfind(1,linei)        = 1
         f_beg(linei,iorder,1) = lpix
         f_end(linei,iorder,1) = upix
         sf_flag(iorder,linei) = .true. 


c     leap frog the feature for continued searching; we increment to the
c     pixel following the upper pixel of the feature and continue
c     scanning the spectrum

         ipix = upix + 1

c     no detection, increment the center pixel and continue scanning

       ELSE  

         ipix = ipix + 1

       END IF

 11   CONTINUE

c     if the ion is HI, then there can be very extended wings that can
c     be parzed into multiple regions by noise.  for Lyalpha, we check
c     if these regions are separated by more than a single resolution
c     element.  if so, we keep them as are, if not, we combine them

C     NOT IMPLEMENTED
C      IF (int(wave0(1,1)).eq.1215) call telescope

      sysfeat = linei 

c     return

      RETURN
      END


c
c.........................................................................
c

      INTEGER FUNCTION subfeat(iorder,linei)

c     This function returns the number of significant absorption
c     subregions within the master region defined by LINEI and the
c     starting and ending pixels, stored in f_beg and f_end, of each
c     feature.  The method uses the equivalent width spectrum, not the
c     simple flux decrement.  Though it can be used for emission
c     features to, this subroutine is designed to scope out absorption
c     features only

c     A feature us found when one pixel has an EW in its pixel that is
c     Nsig beyond the EWSIG uncertainty in that pixel, i.e., EW/EWSIG <
c     -Nsig.  Then the spectrum is scanned backwards to find the blue
c     edge of the feature, which is defined when the ratio EW/EWSIG = -
c     1.  This is the location at about which the flux decrement becomes
c     consistent with the continuum.  Then the spectrum is set back to
c     the initial pixel of the detection and scanned forward for
c     EW/EWSIG = -1.  The pixel locations of the lower and upper bounds
c     are stored in F_BEG and F_END.  Then the pixel number is set to
c     the upper pixel + 1 of the just found feature and the scanning for
c     another Nsig detection is searched for.

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            ipix,lpix,upix,pix1,npix,subregi
      integer            iorder,linei
      double precision   Nsig,ratio,v1,v2,vb,ve

      include            'sysanal.com'



      subregi    = 0
      subfeat = 0

c     set the significance level 

      Nsig  = N_sigma(iorder)

c     find the beginning and ending pixels for the master region
c     (obtained from the first ion in the list) the values PIX1 and NPIX
c     are unmodified in this routine

      pix1 = 0
      npix = 0

      vb = vel(f_beg(linei,1,1),1)
      ve = vel(f_end(linei,1,1),1)
      DO 09 ipix=1,ndata(iorder)-1
       v1 = vel(ipix,iorder)
       v2 = vel(ipix+1,iorder)
       IF ((vb.ge.v1).AND.(vb.le.v2)) pix1 = ipix
       IF ((ve.ge.v1).AND.(ve.le.v2)) npix = ipix
 09   CONTINUE

c     there are now two possibilities if we failed to grab PIX1 and NPIX

c     condition 1: partial data in this region; find the data edges;
c     first check blue edge, then check red edge

      IF ((pix1.eq.0).AND.(npix.ne.0)) pix1 = 1
      IF ((npix.eq.0).AND.(pix1.ne.0)) npix = ndata(iorder)

c     condition 2: data do not extend over range in both blue and red
c     directions or no data in this region at all; which is it?  set
c     accordingly

      IF ((pix1.eq.0).AND.(npix.eq.0)) then
       IF ((vel(1,iorder).gt.vb).AND.(vel(ndata(iorder),iorder).lt.ve)) 
     &   then
        pix1 = 1
        npix = ndata(iorder)
       END IF 

c     no data?; we are too red

       IF (vel(ndata(iorder),iorder).lt.vb) then
        subfeat = 0
        RETURN
       END IF 

c     no data?; we are too blue

       IF (vel(1,iorder).gt.ve) then
        subfeat = 0
        RETURN
       END IF 


      END IF

c     we scan the spectrum over the projected master region and find the
c     subregions for this ion transitions; ipix is the pixel meeting the
c     detection criterion; then we scan blueward and then redward
c     looking for the extremities of the feature; we use a while loop
c     because we need to leap frog the feature once we have it defined,
c     thus we do not increment ipix uniformily

      ipix     = pix1

      DO 11 WHILE (ipix.lt.npix)

       ratio = ewpix(ipix,iorder)/ewsigpix(ipix,iorder)
       
c     enter scan feature block if we have a significant equivalent
c     width; increment the feature index and then scan for the feature
c     extremities

       IF (ratio.le.-Nsig) then 

         subregi = subregi + 1
         lpix    = ipix
         upix    = ipix

c     search blueward for the beginning of the feature

         DO 13 WHILE ((ratio.lt.-1.0).AND.(lpix.gt.1))
           lpix  = lpix - 1 
           ratio = ewpix(lpix,iorder)/ewsigpix(lpix,iorder)  
 13      CONTINUE

c     reset the ratio and pixel to the scanning start point

         ratio = ewpix(ipix,iorder)/ewsigpix(ipix,iorder) 

c     search redward for the ending of the feature

         DO 14 WHILE ((ratio.lt.-1.0).AND.(upix.lt.npix))
           upix  = upix + 1 
           ratio = ewpix(upix,iorder)/ewsigpix(upix,iorder)  
 14      CONTINUE

c     front edge effects; warning

C         IF (lpix.eq.pix1) then 
C           WRITE(6,*) ' WARNING: subfeature',subregi,' artificially ',
C     &                '-terminated start of region ',linei,
C     &                '-for ',order(iorder)(1:20)
C         END IF

c     back edge effects; warning

C         IF (upix.eq.npix) then 
C           WRITE(6,*) ' WARNING: subfeature',subregi,' artificially ',
C     &                '-terminated end of region ',linei,
C     &                '-for ',order(iorder)(1:20)
C         END IF

c     got one- store the pixel indices and set the SL_FLAG high

         f_beg(linei,iorder,subregi) = lpix
         f_end(linei,iorder,subregi) = upix
         sf_flag(iorder,linei) = .true.
         
c     leap frog the feature for continued searching; we increment to the
c     pixel following the upper pixel of the feature and continue
c     scanning the spectrum

         ipix = upix + 1

c     no detection, increment the center pixel and continue scanning

       ELSE  

         ipix = ipix + 1

       END IF

 11   CONTINUE

c     if we have no detection for this ion transition withinin this
c     master region then we set the region size to that of the master
c     region and set the SF_FLAG low; this will allow tracking of
c     limits for further computations


      IF (subregi.eq.0) then
       subregi = 1
       f_beg(linei,iorder,subregi) = pix1
       f_end(linei,iorder,subregi) = npix
       sf_flag(iorder,linei) = .false.
      END IF

      subfeat = subregi

c     return

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE spurious(k,action,flag)

c     eliminate spurious master regions; this routine must be called
c     before the calculations on the data or there will be trouble


c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i,j,k,action
      logical            flag
      double precision   total,knt

      include            'sysanal.com'


c     K is the line (master region number) being examined

      flag   = .false.
      total  = 0         ! = total master regions
      knt    = 0         ! = count of significant regions in higher ions

c     skip the first ion transition; for it does have features!

       DO 11 i=2,norders
        DO 13 j=1,nfind(i,k)

         total = total + 1
         IF (.not.sf_flag(i,k)) knt = knt + 1

 13     CONTINUE
 
 11    CONTINUE

c     if flag is high and action = 1, then delete; if action is low,
c     return flag high

       IF (total.eq.knt) flag = .true.
       IF ((action.eq.1).AND.(flag)) CALL delregion(k)


      RETURN
      END


c

c.........................................................................
c

      SUBROUTINE delregion(linei)

c     telescope the arrays for the removal of a master region

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      include            'sysanal.h'

      integer            i,j,k,linei

      include            'sysanal.com'


       DO 09 i=1,norders

        DO 11 k=linei+1,nlines

         DO 13 j=1,nfind(i,k)

          f_beg(k-1,i,j) = f_beg(k,i,j)
          f_end(k-1,i,j) = f_end(k,i,j)
        
 13      CONTINUE

        nfind(i,k-1)   = nfind(i,k)
        sf_flag(i,k-1) = sf_flag(i,k)

 11     CONTINUE

 09    CONTINUE

       nlines = nlines - 1


       RETURN
       END


