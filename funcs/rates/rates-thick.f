c
c.........................................................................
c

      SUBROUTINE mkfgrid(ilev,kdo,jdo,Zmet,toler)

c     treat self-shielding; 

c     if we enter this routine, then some cells in level ILEV are
c     optically thick; we construct a grid of the mean ionization
c     fractions for the target ion (kdo,jdo) as a function of nH and T.
c     The mean ionization fraction is a weighted mean ionization
c     fraction, where the weights are the fractional pathlengh for each
c     sublayer of the a gas structure with physical depth for the cell
c     size of level ILEV.

c     this routine computes the ionization structure of the model cells
c     and compute the mean ionization fractions and number densities; by
c     model cell, we mean a cell of physical depth given by the size of
c     the cells for level ILEV at the nH,T grid point

c     self shielding is handled by attenuating the ionizing SED as it
c     penetrates into the model cell by dividing the model cell into
c     sublayers (indexed as ISUB); the attenutation of the ionizing SED
c     is treated as a function of optical depth at the HI, HeI, and HeII
c     ionization edges; since the optical depth depends upon the density
c     of the H and He ions, which depsends upon theie ionization
c     fractions, we must solve the ionization balance in each sublayer

c     note, since the grid is rectangular, not all model cells will be
c     optically thick; so some nH,T grid points will ionization
c     fractions for optically thin cells

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c


      implicit none
      include           'rates.h'
      include           'getcube.h'
      logical           opticallythick
      integer           ilev,isub
      integer           in,iT,m,kdo(Imax),jdo(Imax)
      external          meanfreepath
      double precision  zbrent,meanfreepath,dStoler
      parameter         (dStoler=1.0d-5) ! dS accuracy [zbrent]
      double precision  fm(Imax)
      double precision  Lcelli,toler,Zmet,T
      double precision  depth,fdepth,dS,fdS
      double precision  tausum
      character*80      flabel(Imax)
      include           'rates.com'
      include           'com-photo.com'
      include           'com-Rkjs.com'
      include           'getcube.com'


c     create labels for communication of progress (in tabular form) to
c     screen and to the runlog file

      DO m=1,noutfiles
       IF (kdo(m).eq.0) then
        flabel(m) = 'X'   ! if electron density make dummy label
       ELSE
        flabel(m) = 'f'//ionID(kdo(m),jdo(m))
       ENDIF
      ENDDO

c     we use the half-cell depth because we assume symmetry, meaning
c     that the incident ionizing spectrum is incident on two opposing
c     faces of the "plane parallel" cell

      Lcelli  = 0.50d0*Lcell(ilev) ! in pc

c     loop over the density grid

      DO 01 in=1,Nrn

c     stamp header for progress table

      WRITE(6,599) 'inH','iT','Nsub','taumax','lognH','logT',
     &             (flabel(m),m=1,noutfiles)
      WRITE(4,599) 'inH','iT','Nsub','taumax','lognH','logT',
     &             (flabel(m),m=1,noutfiles)

c     set nH and determine the abundances

       nH = 10.0d0**nR(in)
       CALL initabund(Zmet)      

c     loop over the temperature, and initialize all conditions to the
c     optically thin condition for the first layer (the cell face)

       DO 02 iT=1,NRT

        T = 10.0d0**TR(iT)

        isub    = 1
        taumin  = tau(isub)    ! reset optical depth to the cell face
        tausum  = taumin       ! being summation for total optical depth
        CALL mkSED             ! reset SED to the cell face
        CALL pop_photoRs       ! repopulate photo/Auger rates to cell face
        CALL alphas_T(T)       ! compute collisional rate coeffs

c     solve the ionization balance at the cell face; determine whether
c     the model cell (nH,T grid point) is optically thick; if not, store
c     the optically thin result, and jump to the next point on the nH,T
c     grid

        CALL ionbalance(toler) 

        IF (.not.opticallythick(ilev)) then
         DO m=1,noutfiles
          fm(m) = fion(kdo(m),jdo(m))
         ENDDO
         tausum = 0.0d0
         fdepth = 0.0d0
         GOTO 03 ! jump to next point on grid
       ENDIF

c     so, the model cell is optically thick, we boldly go...

c     MODEL CELL FACE RESULTS 

c     compute the model cell face physical depth (dS), effectively the
c     mean free path of the photons for the dominant ionization edge,
c     and the fraction into the half-cell it extends (fdepth)

        dS      = zbrent(meanfreepath,0.0d0,Lcelli,dStoler)

c     begin the bookeeping of the model cell penetration, store the
c     current depth, the fraction depth of the current layer (the model
c     cell face), fdS, and the cuurent fractional depth into the model
c     cell, fdepth

        depth   = dS
        fdS     = dS/Lcelli
        fdepth  = fdS

c     for the target ions, store the first term in the temporary
c     weighted summation, fm(m), m=1,noutfiles, used to obtain the
c     weighted mean ionization fractions for target species kdo,jdo; we
c     weight each term in the sum by fractional pathlength into the
c     model cell, fdS, so that the weighted average ionization fraction
c     preserves the column density of a line of sight through the cell;
c     when we have modelled the entire cell, the value of fdepth, which
c     is the sum of fdS, is 1

        DO m=1,noutfiles
         fm(m) = fion(kdo(m),jdo(m))*fdS
        ENDDO
 
c     MODEL CELL PENETRATION: BEGIN ITERATIVE LOOP FOR SELF-SHIELDING
c     FOR SUBSEQUENT SUBCELLS: we penetrate into the cell one subcell at
c     a time until we reach the physical half-depth of the cell or the
c     maximum optical depth

        DO WHILE ((isub.lt.NRtaumax).AND.(fdepth.lt.1.0))

         isub = isub + 1   ! increment the subcell index counter

c     attenuate the ionizing SED for having penetrated into the
c     *previous* subcell a distance dS; compute photionization and Auger
c     rates and solve ionization balance for this subcell

         CALL attenuateSED(dS)
         CALL alphas_nonT           ! this is the time consuming step
         CALL ionbalance(toler)

c     determine this subcell's physical depth, dS and what fraction of
c     the half cell size, fdS

         dS     = zbrent(meanfreepath,0.0d0,Lcelli,dStoler)
         fdS    = dS/Lcelli

c     if this is the last subcell then fdepth could exceed unity, so
c     correct for that; or, if it is not the last subcell, increment
c     fdepth and depth

         IF ((fdepth+fdS).gt.1.0d0) then
          dS     = Lcelli - depth
          fdS    = dS/Lcelli
          fdepth = 1.0d0
          depth  = depth + dS
         ELSE
          fdepth = fdepth + fdS
          depth  = depth + dS
         END IF

c     increment the next term in the weighted sum used to obtain the
c     weighted mean ionization fractions for the target ions

         DO m=1,noutfiles
          fm(m) = fm(m) + fion(kdo(m),jdo(m))*fdS
         ENDDO

c     set the optical depth increment taumin (global variable) for the
c     attenuation of this layer as it applies to the next layer;
c     increment the total optical depth (used only for reporting)

         taumin = tau(isub)
         tausum = tausum + tau(isub)

        ENDDO  ! end WHILE loop- step to the next sub cell

c     SELF-SHIELDING TREATMENT COMPLETE OR OPTICAL DEPTH MAXIMUM REACHED

c     TRAPPING HIGHER OPTICAL DEPTH GRID POINTS: if this model cell is
c     super optically thick, we stopped at isub=NRtaumax, so we project
c     the next layer's ionization fraction at this maximum depth to
c     treat the remainder of the interior of the cell; this also insures
c     that the sum of fdS equals unity so that the weighting is done
c     proper; NOTE: this assumes that the remainder of the model cell is
c     unchanging in its ionization conditions, which is a pretty good
c     assumption (from testing) and the fact that at this very high
c     optical depth, the photoionization rates are not changing
c     substantially because they have become negligible compared to
c     collisional processes, though for very low temperatures this
c     approximation can break down because the collisional rate
c     coefficients can vanish so that the cell is dominated by
c     photoionization

        IF (isub.eq.NRtaumax) then ! project to back edge of model cell
         CALL attenuateSED(dS)
         CALL alphas_nonT
         CALL ionbalance(toler)
         dS     = Lcelli - depth
         fdS    = dS/Lcelli
         fdepth = 1.0d0
         DO m=1,noutfiles
          fm(m) = fm(m) + fion(kdo(m),jdo(m))*fdS
         ENDDO
        END IF

c     we are done for this grid point; store the weighted mean
c     ionization fractions at the nH,T grid point for the target ions;
c     if the grid location is optically thin, then we jumped here from
c     above the WHILE loop

 03     DO m=1,noutfiles
         fgrid(m,in,iT) = fm(m)
        ENDDO

c     report the structure of this cell and the weighted ionization
c     fraction for the target ions

        WRITE(6,600) in,iT,isub,tausum,nR(in),TR(iT),
     &               (fm(m),m=1,noutfiles)
        WRITE(4,600) in,iT,isub,tausum,nR(in),TR(iT),
     &               (fm(m),m=1,noutfiles)

 02    CONTINUE ! next T in the grid
 01   CONTINUE  ! next nH in the grid

      RETURN

 599  FORMAT(1x,3a5,3a12,4x,20a10) 
 600  FORMAT(1x,3i5,3f12.3,1x,1p20e10.2) 

      END

c
c.........................................................................
c

      SUBROUTINE selfshield(T,kdo,jdo)

c     employs the ionization fraction grid to obtain the self-shielding
c     mean ionization condition for the current cell

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c


      implicit none
      include           'rates.h'
      include           'getcube.h'
      integer           in,iT,k,m,kdo(Imax),jdo(Imax)
      integer           knlo,knhi,kTlo,kThi
      double precision  T,lgT,lgnH,hn,hT
      double precision  a,b,c,d
      include           'rates.com'
      include           'com-Rkjs.com'
      include           'getcube.com'


c     convert nH and T into log coordinates for the grid matching

      lgnH = log10(nH)
      lgT  = log10(T)

c     bracket the density location on the grid

      knlo = 1
      knhi = NRn
01    IF (knhi-knlo.gt.1) THEN
        k = (knhi+knlo)/2 
        IF (nR(k).gt.lgnH) THEN
          knhi = k
        ELSE
          knlo = k
        ENDIF
      GOTO 01
      ENDIF

      hn = nR(knhi) - nR(knlo)
      IF (hn.EQ.0.0) then
       WRITE(6,*) 'ERROR(selfshield): nH was not bracketed'
       WRITE(6,*) 'nH: ',lgnH,knlo,knhi,nR(knhi),nR(knlo)
       STOP
      ENDIF

c     bracket the temperature location on the grid

      kTlo  = 1
      kThi  = NRT
 02   IF (kThi-kTlo.GT.1) THEN
        k = (kThi+kTlo)/2 
        IF (TR(k).GT.lgT) THEN
          kThi = k
        ELSE
          kTlo = k
        ENDIF
      GOTO 02
      ENDIF

      hT = TR(kThi) - TR(kTlo)
      IF (hT.EQ.0.0) then
       WRITE(6,*) 'ERROR(selfshield): T was not bracketed'
       WRITE(6,*) 'T : ',lgT,kTlo,kThi,TR(kThi),TR(kTlo)
       STOP
      ENDIF

c     apply bilinear interpolation on the grid to obtain that target ion
c     ionization fractions

      DO m=1,noutfiles
       IF (kdo(m).ne.0) then ! do not do the electron density
        a = fgrid(m,knlo,kTlo)*(nR(knhi)-lgnH)*(TR(kThi)-lgT)
        b = fgrid(m,knhi,kTlo)*(lgnH-nR(knlo))*(TR(kThi)-lgT)
        c = fgrid(m,knlo,kThi)*(nR(knhi)-lgnH)*(lgT-TR(kTlo))
        d = fgrid(m,knhi,kThi)*(lgnH-nR(knlo))*(lgT-TR(kTlo))
        fion(kdo(m),jdo(m)) = (a+b+c+d)/(hn*hT) 
        nkj(kdo(m),jdo(m))  = fion(kdo(m),jdo(m))*nk(kdo(m))
      ENDIF
      ENDDO
  
      RETURN
      END


c
c.........................................................................
c

      DOUBLE PRECISION FUNCTION meanfreepath(dS)

c     this function is called by zebrent, which finds the incremental
c     path length dS for optical depth taumin (global variable)
      
c     function meanfreepath returns dS in units of parsecs

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c


      implicit none
      include           'rates.h'
      include           'getcube.h'
      double precision  dS,chi11,chi21,chi22,chimax
      include           'rates.com'
      include           'com-photo.com'
      include           'com-Rkjs.com'
      include           'getcube.com'

c     compute the opacities at the neutral hydrogen edge (11), the
c     neutral helium edge (21), and the singly ionized helium edge (22);
c     opacity units are inverse centimeters and are converted to inverse
c     parsecs so that dS is returned in parsecs

      chi11 = (nkj(1,1)*sig11) * 3.08567758d+18 
      chi21 = (nkj(2,1)*sig21) * 3.08567758d+18 
      chi22 = (nkj(2,2)*sig22) * 3.08567758d+18 

c     find the dominant ionization edge from the largest opacity

      chimax = max(chi11,chi21)
      chimax = max(chimax,chi22)

c     zero the function such that the value of dS yields optical depth
c     taumin=chimax*dS

      meanfreepath = taumin - (chimax)*dS 

      RETURN
      END

c
c.........................................................................
c

      LOGICAL FUNCTION opticallythick(ilev)

c     global boolean OPTICALLYTHICK is returned high if the mean free
c     path of any of ionizing photons are less than half the cell
c     thickness, and is returned low otherwise; if high, self-shielding
c     will be treated 

c     for cell level ILEV this routine computes the mean free path of
c     the photons at the ionization edges for unity optical depth for
c     neutral hydrogen, neutral helium, or singly ionized helium; then
c     the minimum mean free path from any one of these edges is
c     determined and applied to check the criterion for optical
c     thickness

c     this routine is called one cell level at a time, all cells in the
c     level are checked (only once per cell) in order to determine the
c     nH and T ranges for which self-shielding treatment is required for
c     this cell level

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c


      implicit none
      include           'rates.h'
      include           'getcube.h'
      integer           ilev
      double precision  mfp11,mfp21,mfp22,mfp
      include           'rates.com'
      include           'com-photo.com'
      include           'com-Rkjs.com'
      include           'getcube.com'

c     initialize the function 

      opticallythick = .false.

c     compute the mean free path for unity optical depth at neutral
c     hydrogen edge (11), at the neutral helium edge (21), and at the
c     singly ionized helium edge (22); the mean free paths are converted
c     to parsecs from centimeters

      mfp11 = 1.0d0/(nkj(1,1)*sig11)/3.08567758d+18 
      mfp21 = 1.0d0/(nkj(2,1)*sig21)/3.08567758d+18
      mfp22 = 1.0d0/(nkj(2,2)*sig22)/3.08567758d+18

c     determine which mean free path yielding the the smallest
c     penetration into the cell

      mfp = min(mfp11,mfp21)
      mfp = min(mfp22,mfp)

c     using the half-thickness of the cell for this level, check whether
c     this cell is optically thick; note that we assume that the cell
c     ionization structure is symmetric, or basically that the cell is
c     illuminated on both sides; we will solve the radiative transfer
c     through half the cell and reflect it about the midplane, thus the
c     criterion for ionization structure is that the mean free path is
c     less than the half-depth of the cell

      IF (mfp.lt.0.50d0*Lcell(ilev)) opticallythick = .true.

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE attenuateSED(dS)

c     this routine takes the current SED and attenuates for the
c     incremental pathlength dS into the cell

c     dS is input to the routine in units of parsecs and must be
c     converted to centimeters as necessary

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c


      implicit none
      include           'rates.h'
      include           'getcube.h'
      integer           ie
      double precision  dS,phfit2
      double precision  E,logJE0,y
      double precision  sig11E,sig21E,sig22E
      double precision  tau11E,tau21E,tau22E
      include           'rates.com'
      include           'com-photo.com'
      include           'com-Rkjs.com'


c     attenuate the spectrum for hydrogen and helium ionization; the
c     attenuation is performed in log10 space
 
      DO 31 ie=1,NEsed
       tau11E = 0.0d0
       tau21E = 0.0d0
       tau22E = 0.0d0
       IF (logEsed(ie).ge.log10(E11)) then 
        E      = 10.0d0**logEsed(ie)
        logJE0 = logJEsed(ie)
        sig11E = 1.0d-18*phfit2(1,1,1,E)
        sig21E = 1.0d-18*phfit2(2,2,1,E)
        sig22E = 1.0d-18*phfit2(2,1,1,E)
        tau11E = nkj(1,1)*sig11E*dS*3.08567758d+18 
        tau21E = nkj(2,1)*sig21E*dS*3.08567758d+18 
        tau22E = nkj(2,2)*sig22E*dS*3.08567758d+18 
        y      = 0.434294481d0*(tau11E+tau21E+tau22E)
        logJEsed(ie) = logJE0 - y
       END IF
31    CONTINUE

c     obtain the global 2nd derivatives in the energy direction of
c     logJEuvb for the integration required to obtain the
c     photoionization and Auger rates (performed elsewhere)

      CALL SEDderivs

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE popedges

c     this subroutine populates the photoionization cross sections at
c     the threshold energies for the ionization edges of neutral
c     hydrogen (11), neutral helium (21), and single ionized helium
c     (22); they are stored globally for checking the mean free paths of
c     the ionizing SED 

c     the cross sections stored as cm^2; since there is only a single
c     shell for these ions, there is no summing over shells

c     function thresholdE returns the threshold energy E for ionization

c     function phfit2 returns the cross section in unit of 10^18 cm^2 at
c     photon energy E (i.e., threshold energies E11, E21, E22)
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      double precision  thresholdE,phfit2
      include           'com-photo.com'

c     get the cross section at neutral hydrogen edge

      E11   = thresholdE(1,1,1) 
      sig11 = 1.0d-18 * phfit2(1,1,1,E11)

c     get the cross section at neutral helium edge

      E21   = thresholdE(2,2,1)
      sig21 = 1.0d-18 * phfit2(2,2,1,E21)

c     get the cross section at ionized helium edge

      E22   = thresholdE(2,1,1)
      sig22 = 1.0d-18 * phfit2(2,1,1,E22)

      RETURN
      END

