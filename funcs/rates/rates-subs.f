c.........................................................................
c

      SUBROUTINE configall(toler,GZinfile,z)

c     this routine configures the run, including 
c     - read files containing book keeping data
c     - initialize (zero) array/matrices and setup book keeping array
c     - determine the abundance fractions
c     - read in fitting parameters for the rate coefficients
c     - compute the ionizing spectrum (SED)
c     - set physical sizes of arrays and set array pointers

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      double precision  toler,z
      character*80      GZinfile  ! gas+metallicity file
      include           'rates.com'


c     ........................
c     BOOK KEEPING TABLE READS
c     ........................

c     read a few tables of information for user interfacing ease
c     - the list of included species, indexed by k
c     - ion IDs (i.e., MgII, etc) which are indexed by k,j
c     - shell IDs (as defined by Verner) indexed by s;
c     - electron confgurations for iso-electronic sequences (1s1, etc)

      CALL inputdek(toler,GZinfile,z)
      CALL readionIDs
      CALL readshellIDs
      CALL readeconfigs

c     Note: global variable Nspecies returned by routine readionIDs


c     .................
c     INITIALIZE ARRAYS
c     .................

c     initialize (null) all fitting parameters (routine zeroall);
c     initialize book keeping arrays for charge exchange recombinations
c     with He (routine Hesetup)

      CALL zeroall
      CALL Hesetup


      RETURN

      END


c.........................................................................
c

      SUBROUTINE configIonPhysics(z)

c     communicate the physics employed and load in the needed data

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      logical           doCOLL
      double precision  z
      include           'rates.com'
      include           'com-modes.com'

c     check if any collisional processes are included 

      doCOLL = .false.
      IF (doLDR.OR.doHDR.OR.doCDI.OR.doCEA.OR.doCT) doCOLL = .true.

c     begin communication of the included physics

      CALL notify(' ',-1)
      CALL notify(' IONIZATION PHYSICS',-1)

c     sanity check to avoid crashing

      IF ((.not.doPH).AND.(.not.doCOLL)) then
       CALL notify(' *********************************',-1)
       CALL notify(' *   FATAL CONFIGURATION!        *',-1)
       CALL notify(' *   Both photo and collisional  *',-1)
       CALL notify(' *   ionization are turned OFF.  *',-1)
       CALL notify(' *   You must include some form  *',-1)
       CALL notify(' *   of ionization process!!     *',-1)
       CALL notify(' *   Reconfigure the inputs in   *',-1)
       CALL notify(' *   the rates.inp file.         *',-1)
       CALL notify(' *********************************',-1)
       STOP
      ENDIF

      IF ((.not.doPH).AND.(.not.uvbflag).AND.(.not.sb99flag)) then
       CALL notify(' *********************************',-1)
       CALL notify(' *   FATAL CONFIGURATION!        *',-1)
       CALL notify(' *   Photoionization is ON, but  *',-1)
       CALL notify(' *   the ionizing SED is not     *',-1)
       CALL notify(' *   specified.                  *',-1)
       CALL notify(' *   Reconfigure the inputs in   *',-1)
       CALL notify(' *   the rates.inp file.         *',-1)
       CALL notify(' *********************************',-1)
       STOP
      ENDIF

c     photoionization?

      IF (doPH) then

       CALL notify('  * photoionization is ON',-1)

c     self shielding?
       IF (doSLFSH) then                                     
        CALL notify('  * self-sheilding is ON',-1)
       ELSE
        CALL notify('  * self-shielding is OFF',-1)
       ENDIF

c     SED usage
       IF ((uvbflag).AND.(.not.sb99flag)) then                
        CALL notify('   + photoionization is UVB only',-1)
       ENDIF        
       IF ((.not.uvbflag).AND.(sb99flag)) then
        CALL notify('   + photoionization is stars only',-1)
       ENDIF        
       IF ((uvbflag).AND.(sb99flag)) then
        CALL notify('   + photoionization is UVB + stars',-1)
       ENDIF        

c     loading SED libraries
       IF (uvbflag) then
        CALL notify('     loading UVB SED library',-1)
        CALL readUVBspec
        CALL mkUVBofz(z)
       ENDIF
       IF (sb99flag) then
        CALL notify('     loading Starburst 99 libraries',-1)
        CALL readSB99spec  
       ENDIF

c     loading photo cross sections
       CALL notify('     loading photo cross-section data',-1)
       CALL popedges
       CALL photoxsecs

      ELSE  ! not doing photoionization

       CALL notify('  * photoionization is OFF',-1)

      ENDIF 

c     Auger ionization

      IF (doAUG) then
       CALL notify('  * Auger ionization is ON',-1)
       CALL notify('     loading Auger yields',-1)
       CALL readAugerKM93
      ELSE
       CALL notify('  * Auger ionization is OFF',-1)
      ENDIF

c     collisional ionization

      IF (doCDI.OR.doCEA.OR.doCT) then
       CALL notify('  * collisional ionization is ON:',-1)
       IF (doCDI) CALL notify('   + direct collisional ionization',-1)
       IF (doCEA) CALL notify('   + excitation autoionization',-1) 
       IF (doCT)  CALL notify('   + charge-exchange ionization',-1) 
       IF (doCDI.OR.doCEA) then 
        CALL notify('     loading collisional cross-section data',-1)
        CALL readCdiA85tab 
       ENDIF
      ELSE
       CALL notify('  * collisional ionization is OFF:',-1)
      ENDIF 

c     recombination (radiative is ALWAYS on)

       CALL notify('  * two-body radiative recombination is ON:',-1)
       IF (doLDR) then
        CALL notify('  * low-T dielectronic recombination is ON',-1)
       ELSE
        CALL notify('  * low-T dielectronic recombination is OFF',-1) 
       ENDIF
       IF (doHDR) then
        CALL notify('  * high-T dielectronic recombination is ON',-1) 
       ELSE
        CALL notify('  * high-T dielectronic recombination is OFF',-1) 
       ENDIF
       IF (doCT) then
        CALL notify('  * charge-exchange recombination is ON',-1) 
       ELSE
        CALL notify('  * charge-exchange recombination is OFF',-1) 
       ENDIF
       CALL notify('     loading recombination data',-1)
       CALL readPLtab
       IF (doLDR) CALL readLDRtab
       IF (doHDR) CALL readHDRtab


      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE inputdek(toler,GZinfile,z)

c     this routine reads the file 'rates.inp' 

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      integer           i,j,k,nel
      double precision  z,toler,D1MACH,afac
      character*5       flag
      character*80      header
      character*80      GZinfile  ! gas+metallicity file
      include           'rates.com'
      include           'com-modes.com'


c     open the 'rates.inp' file and set the inclusion flags for the
c     species

      OPEN(unit=1,file='rates.inp',status='old')

      CALL notify(' GATHERING INPUT CONFIGURATION',-1)

c     strip off the file master header

      READ(1,*) header

c     ............
c     RUN MODE
c     ............

      CALL notify('  obtaining box info...',-1)

c     strip off header

      READ(1,*) header
      READ(1,*) header

      READ(1,*) GZinfile
      READ(1,*) afac
      READ(1,*) noutfiles

      z = 1.0d0/afac - 1.0d0
      IF (z.lt.0.0d0) z = 0.0d0

c     the run mode (integer)

c     ............
c     RUN MODE
c     ............

      CALL notify('  obtaining runtime mode...',-1)

c     strip off header

      READ(1,*) header
      READ(1,*) header

c     UVB inclusion boolean

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'yes'.OR.flag.eq.'YES') uvbflag = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   uvbflag = .false.

c     SB99 inclusion boolean

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'yes'.OR.flag.eq.'YES') sb99flag = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   sb99flag = .false.

c     ........
c     LSF MODE
c     ........

      CALL notify('  obtaining solution mode...',-1)

c     strip off header

      READ(1,*) header
      READ(1,*) header

c     obtain the self-shielding flag

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'no') doSLFSH = .false.
      IF (flag.eq.'yes'.OR.flag.eq.'YES') doSLFSH = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   doSLFSH = .false.

c     obtain the desired tolerance to the fit

      READ(1,*) toler
      IF (toler.eq.0.) toler = SQRT(D1MACH(4))


c     ............
c     PHYSICS MODE
c     ............
c     note: photo-recombination must always included 

      CALL notify('  obtaining ionization physics modes...',-1)

c     strip off header

      READ(1,*) header
      READ(1,*) header

c     (PH) photoionzation boolean

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'yes'.OR.flag.eq.'YES') doPH = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   doPH = .false.

c     (LDR) low temperature dielectronic recombination boolean

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'yes'.OR.flag.eq.'YES') doLDR = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   doLDR = .false.

c     (HDR) high temperature dielectronic recombination boolean

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'yes'.OR.flag.eq.'YES') doHDR = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   doHDR = .false.

c     (CDI) direct collisional ionization boolean
c     if PH is turned off, we turn CDI on

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'yes'.OR.flag.eq.'YES') doCDI = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   doCDI = .false.
      IF (.not.doPH) doCDI = .true.

c     (CEA) collisional excitation-autoionization boolean
c     if CDI is turned off, we turn off CEA

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'yes'.OR.flag.eq.'YES') doCEA = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   doCEA = .false.
      IF (.not.doCDI) doCEA = .false.

c     (AUG) auger ionzation boolean
c     this is photoionization effect; if PH is turned off then
c     this gets turned off too

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'yes'.OR.flag.eq.'YES') doAUG = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   doAUG = .false.
      IF (.not.doPH) doAUG = .false.

c     (CT) charge transfer boolean

      READ(1,'(a3)') flag
      IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
      IF (flag.eq.'yes'.OR.flag.eq.'YES') doCT = .true.
      IF (flag.eq.'no'.OR.flag.eq.'NO')   doCT = .false.



c     .................
c     METALS TO INCLUDE
c     .................

      CALL notify('  obtaining the chemical mixture...',-1)

c     we always do hydrogen and helium

      specflag(1) = .true.
      specflag(2) = .true.      

c     strip the header

      READ(1,*) header
      READ(1,*) header

c     read in the metal inclusion booleans

      DO 11 k=3,Imax
       specflag(k) = .false.
       READ(1,'(a3)',END=12) flag
       IF ((flag.ne.'yes').AND.(flag.ne.'YES').AND.
     &     (flag.ne.'no').AND.(flag.ne.'NO')) GOTO 14
       IF (flag.eq.'yes'.OR.flag.eq.'YES') specflag(k) = .true.
 11   CONTINUE

      CLOSE(unit=1)
      RETURN

 12   CLOSE(unit=1)
      WRITE(6,*) 'ERROR(inputdek): end of file reached'
      STOP

 14   CLOSE(unit=1)
      WRITE(6,*) 'ERROR(inputdek): bad entry in rates.inp file'      
      WRITE(6,*) 'specify "yes" or "no" for boolean flags.'      
      STOP

      END

c
c.........................................................................
c

      SUBROUTINE initabund(Zmet)

c     populates the species mass and abundance fractions and number
c     densities

c     abundance data taken from Drain 2011, "Physics of the Interstellar
c     and Intergalactic Medium", p8 Table 1.4; all data lifted from
c     original work of Asplund et al (2009,ARA&A,47,481)

c     atomic masses (in amu) taken from NIST

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           k,kk
      double precision  ZoXtot,ZoXinc,Zmet,asum,xsum
      double precision  xkxH(Imax),amu(Imax),xr
      double precision  xksol(Imax)
      include           'rates.com'


c     populate the primordial solar abundance mass fraction ratios
c     relative to H (Drain 2011)

c                xk/xH    Elem  Ref

      xkxH(1)  = 1.0d0    ! H
      xkxH(2)  = 3.82d-01 ! He   photospheric
      xkxH(3)  = 1.38d-08 ! Li   meteoric
      xkxH(4)  = 1.97d-10 ! Be   met
      xkxH(5)  = 7.31d-09 ! B    met
      xkxH(6)  = 3.54d-03 ! C    photo
      xkxH(7)  = 1.04d-03 ! N    photo
      xkxH(8)  = 8.59d-03 ! O    photo
      xkxH(9)  = 5.48d-07 ! F    met
      xkxH(10) = 1.88d-03 ! Ne   photo
      xkxH(11) = 4.96d-05 ! Na   met
      xkxH(12) = 1.06d-03 ! Mg   photo
      xkxH(13) = 8.85d-05 ! Al   met
      xkxH(14) = 9.07d-04 ! Si   photo
      xkxH(15) = 1.00d-05 ! P    photo
      xkxH(16) = 4.63d-04 ! S    photo
      xkxH(17) = 6.60d-06 ! Cl   met
      xkxH(18) = 1.10d-04 ! Ar   photo
      xkxH(19) = 5.15d-06 ! K    met
      xkxH(20) = 8.57d-05 ! Ca   met
      xkxH(21) = 5.53d-08 ! Sc   met
      xkxH(22) = 4.27d-06 ! Ti   met
      xkxH(23) = 5.09d-07 ! V    met
      xkxH(24) = 2.49d-05 ! Cr   met
      xkxH(25) = 1.82d-05 ! Mn   met
      xkxH(26) = 1.94d-03 ! Fe   photo
      xkxH(27) = 4.97d-06 ! Co   met
      xkxH(28) = 1.02d-04 ! Ni   met
      xkxH(29) = 1.24d-06 ! Cu   met
      xkxH(30) = 3.06d-06 ! Zn   met


c     atomic mass in amu [NIST]

      amu(1)  =  1.0794d0    ! H
      amu(2)  =  4.002602d0  ! He
      amu(3)  =  6.941d0     ! Li
      amu(4)  =  9.012182d0  ! Be
      amu(5)  = 10.811d0     ! B
      amu(6)  = 12.0107d0    ! C
      amu(7)  = 14.0067d0    ! N
      amu(8)  = 15.9994d0    ! O
      amu(9)  = 18.9984032d0 ! F
      amu(10) = 20.1719d0    ! Ne
      amu(11) = 22.989770d0  ! Na
      amu(12) = 24.3050d0    ! Mg
      amu(13) = 26.981538d0  ! Al
      amu(14) = 28.0855d0    ! Si
      amu(15) = 30.973761d0  ! P
      amu(16) = 32.065d0     ! S
      amu(17) = 35.453d0     ! Cl
      amu(18) = 39.948d0     ! Ar
      amu(19) = 39.0938d0    ! K
      amu(20) = 40.078d0     ! Ca
      amu(21) = 44.955910d0  ! Sc
      amu(22) = 47.867d0     ! Ti
      amu(23) = 50.9415d0    ! V
      amu(24) = 51.9961d0    ! Cr
      amu(25) = 54.938049d0  ! Mn
      amu(26) = 55.845d0     ! Fe
      amu(27) = 58.933200d0  ! Co
      amu(28) = 58.6934d0    ! Ni
      amu(29) = 63.546d0     ! Cu
      amu(30) = 65.409d0     ! Zn

c     for the full mixture up to and including k=30 (Zn), rescale the
c     relative mass fractions from Drain (2011) to absolute mass
c     fractions for the primordial solar abundance pattern

      xsum = 0.0d0
      DO 101 k=1,Imax
       xsum = xsum + xkxh(k)
 101  CONTINUE
      DO 103 k=1,Imax
       xksol(k) = xkxh(k)/xsum  ! solar mass fractions
 103  CONTINUE

c     convert full mixture mass fraction to full mixture abundance
c     fractions

      asum = 0.0d0
      DO 105 k=1,Imax
       asum = asum + xksol(k)/amu(k)
 105  CONTINUE
      DO 107 k=1,Imax
       alphksol(k) = (xksol(k)/amu(k))/asum  ! solar abundance fractions 
 107  CONTINUE

c     obtain the normalizations to scale the solar mass fractions to
c     account for only those species being including in the mixture of
c     the model; ZoXtot is the metal mass fraction for the full mixture;
c     ZoXinc is the metal mass fraction for the included mixture

      ZoXtot = 0.0d0
      DO 01 k=3,Imax
       ZoXtot = ZoXtot + xksol(k)   ! full mixture total
 01   CONTINUE

      ZoXinc = 0.0d0
      DO 03 kk=3,Nspecies
       k = kidx(kk)
       ZoXinc = ZoXinc + xksol(k)   ! included mixture total
 03   CONTINUE

c     rescale the metal mass fractions of the included mixture to the
c     constraint mass fraction, Zmet [*** NOTE COMMENTED OUT LINE ***]

      asum = 0.0d0
      DO 05 kk=3,Nspecies
       k = kidx(kk)
C       xfrac(k) = ((Zmet*ZoXtot)/ZoXinc) * xksol(k) ! If Zmet in solar
       xfrac(k) = (Zmet/ZoXinc) * xksol(k)           ! If Zmet absolute
       asum = asum + xfrac(k)
  05  CONTINUE

c     we must now rescale the hydrogen and helium mass fractions so that
c     the total mass fractions sum to unity (X+Y+Z=1), where Z=ASUM
c     computed in the above loop; we substitute Y=XR*X and solve for X,
c     then for Y.  We adopt the solar mass fraction ratio
c     xr=xkxh(2)/xkxh(1), i.e., (He/H), from the primordial solar
c     mixture (until we get the data for SNII and SNIa mixtures)

      xr       = xkxh(2)/xkxh(1)
      xfrac(1) = (1.0d0-asum)/(1.0d0+xr)  ! solve hydrogen
      xfrac(2) = xr*xfrac(1)              ! solve helium

C     sum check
C      asum = 0.0d0
C      DO kk=1,Nspecies
C       k = kidx(kk)
C       asum = asum + xfrac(k)
C      ENDDO
C      WRITE(6,*) 'sum of mass fractions = ',asum
C      WRITE(6,*) 'hydrogen              = ',xfrac(1)
C      WRITE(6,*) 'helium                = ',xfrac(2)
C      WRITE(6,*) 'metal                 = ',asum-xfrac(1)-xfrac(2)
C      WRITE(6,*) 'should equal          = ',Zmet      
C      STOP

c     now compute final model cloud abundance fractions for the scaled
c     mixture

      asum = 0.0d0
      DO 109 kk=1,Nspecies
       k    = kidx(kk)
       asum = asum + xfrac(k)/amu(k)
 109  CONTINUE
      DO 111 kk=1,Nspecies
       k         = kidx(kk)
       fabund(k) = (xfrac(k)/amu(k))/asum  ! included mixture abund fracs
 111  CONTINUE

c     NOTE: at this point, we have the included mixture mass fractions
c     and abundance ratios; the mass fractions of metals add to
c     Zmet*ZoXtot, the mass fractions of all species add to unity, and
c     the abundance fractionsof all species add to unity [again, I ran
c     tests and have verified this]

c     compute the number densities for each included species in the
c     mixture; they are scaled by the abundance fractions relative to
c     the input hydrogen number density

      nions = 0.0d0
      DO 09 kk=1,Nspecies
       k         = kidx(kk)
       nk(k)     = (fabund(k)/fabund(1)) * nH
       nions     = nions + nk(k)
 09   CONTINUE

      RETURN
      END


c
c.........................................................................
c

      SUBROUTINE zeroall

c     this routine nulls the rate matrix elements (ans several other
c     arrays)

c     it then initializes all arrays for parameters that are read via
c     files, which are the recombinations rate parameters (radiative,
c     loT and hiT dielectronic), the collisional direct ionization rate
c     paramaters, and the Auger yeild probabilities

c     the photoionization rate parameters are initialized in a data
c     block in module phfit-data.blk

c     the charge exchange rate parameters are initialized in a data
c     block in module kingdon-data.blk

c     the excitation-autoionization rates are based upon individual
c     formulae that are directly typed in for each iso-electronic
c     sequence, so there are no arrays to initialize

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           k,j,m,s,i,y
      include           'rates.com'
      include           'com-recomb.com'
      include           'com-coll.com'
      include           'com-auger.com'



c     .......................................................
c     PHOTO RATES, COLLISIONAL RATE COEFFICIENTS, RATE MATRIX
c     .......................................................

c     zero the r

c     k denotes species number
c     j denotes ionzation stage

      DO 01 k=1,Imax                ! loop over species
       DO 02 j=1,Imax+1             ! loop over ionization stages
        R_ph(k,j)          = 0.0d0  ! photo ionization rate
        alpha_rec(k,j)     = 0.0d0  ! total recombination 
        alpha_rec_ph(k,j)  = 0.0d0  ! photo recombination 
        alpha_rec_die(k,j) = 0.0d0  ! dielectronic recombination
        alpha_Cdi(k,j)     = 0.0d0  ! direct collisional ionization
        alpha_Cea(k,j)     = 0.0d0  ! ex-autoionization
        alpha_CTrecH(k,j)  = 0.0d0  ! charge transfer H recombination
        alpha_CTionH(k,j)  = 0.0d0  ! charge transfer H ionization
        alpha_CTrecHe(k,j) = 0.0d0  ! charge transfer He recombination
        R_hat(k,j)         = 0.0d0  ! summed destruction rate of ion k,j
        R_ion(k,j)         = 0.0d0  ! summed ionization rate of ion k,j
        R_rec(k,j)         = 0.0d0  ! summed recombination rate of ion k,j
 02    CONTINUE
 01   CONTINUE

      DO 03 s=1,Nshmax              ! loop over shells
       R_phs(s) = 0.0d0             ! shell photoionization rate
 03   CONTINUE

c     ...................................
c     AUGER RATES (see function AUG_Rate) 
c     ...................................

c     initialize book keeping arrays and zero the fitting parameters for
c     the Auger yield probabilities; note that we initialize the
c     probability W_Auger(k,j,s,1)=1, this is because many shells of an
c     ion do not result in Auger electrons, so that the probability of a
c     single electron ejection from photoionization is unity; in cases
c     where Auger electrons may be ejected, this probability is modified
c     when the probability table is read in (routine readAugerKM93)

c     k = species number
c     j = ionzation stage
c     s = Verner's shell number
c     y = Auger yield (number of electrons ejected)

      DO 31 k=1,Imax               ! loop over species
       DO 32 j=1,Imax              ! loop over ionization stages
        R_Agrout(k,j)  = 0.0d0     ! summed outward Auger rates
        DO 35 i=1,Imax+1           ! loop over final ionization stages
         Q(k,j,i)  = 0.0d0         ! Auger ionization rates k,j >1 electron
 35     CONTINUE
        DO 33 s=1,Nshmax           ! loop over shells (subshells combined)
         W_Auger(k,j,s,1) = 1.0d0  ! non-auger yield probability [DO NOT CHANGE]
         DO 34 y=2,10              ! loop over electron yield
          W_Auger(k,j,s,y) = 0.0d0 ! auger yield probability (Kaastra93)
 34      CONTINUE
 33     CONTINUE
 32    CONTINUE
 31   CONTINUE


c     .....................................
c     RECOMBINATION (see function REC_Rate)
c     .....................................

c     initialize book keeping arrays and zero the fitting parameters for
c     the radiative recombination, low temperature dielectronic
c     recombination, and high temperature dielectronic recombination
c     rate coefficients

c     k denotes species number
c     j denotes ionzation stage
c     i denotes temperature regime (loT dielectronic recomb)

      DO 11 k=1,Imax          ! loop over species
       n1stg(k) = 100         ! loT recomb initial shell (will use min func)
       nNstg(k) = 0           ! loT recomb final shell
       DO 12 j=1,STGmax       ! loop over ionization stages
        nTreg(k,j) = 0        ! loT recomb # of T regions
        A_pl(k,j) = 0.0d0     ! rad recomb fit param A (Arnaud85)
        B_pl(k,j) = 0.0d0     ! rad recomb fit param B (Arnaud85)
        DO 13 i=1,2           ! loop over loT recomb temperature regimes
         a_ldr(k,j,i) = 0.0d0 ! loT recomb fit param a (Nussbaumer96)
         b_ldr(k,j,i) = 0.0d0 ! loT recomb fit param b (Nussbaumer96)
         c_ldr(k,j,i) = 0.0d0 ! loT recomb fit param c (Nussbaumer96)
         d_ldr(k,j,i) = 0.0d0 ! loT recomb fit param d (Nussbaumer96)
         f_ldr(k,j,i) = 0.0d0 ! loT recomb fit param f (Nussbaumer96)
 13     CONTINUE
        A_hdr(k,j)  = 0.0d0   ! hiT recomb fit param A (Arnaud85)
        B_hdr(k,j)  = 0.0d0   ! hiT recomb fit param B (Arnaud85)
        T0_hdr(k,j) = 0.0d0   ! hiT recomb lower T limit (Arnaud85)
        T1_hdr(k,j) = 0.0d0   ! hiT recomb upper T limit (Arnaud85)
 12    CONTINUE
 11   CONTINUE


c     .....................................................
c     DIRECT COLLISIONAL IONIZATION (see function CDI_Rate) 
c     .....................................................

c     initialize book keeping arrays and zero the fitting parameters for
c     the direct collisional ionization rate coefficients

c     k denotes species number
c     j denotes ionzation stage
c     s denotes subshell number (not the same index as Verner's shells)

      DO 21 k=1,Imax          ! loop over species 
       DO 22 j=1,Imax         ! loop over ionization stages    
        shmax(k,j) = 0        ! number of subshells for k,j
        DO 23 s=1,Smax        ! loop over subshells
         IE(k,j,s) = 0.0d0    ! subshell threshold energy I (Arnaud85)
         A(k,j,s)  = 0.0d0    ! di coll ioniz fit param A   (Arnaud85)
         B(k,j,s)  = 0.0d0    ! di coll ioniz fit param B   (Arnaud85)
         C(k,j,s)  = 0.0d0    ! di coll ioniz fit param C   (Arnaud85)
         D(k,j,s)  = 0.0d0    ! di coll ioniz fit param D   (Arnaud85)
 23     CONTINUE
 22    CONTINUE
 21   CONTINUE


c     ............................
c     ABUNDANCE AND GAS QUANTITIES 
c     ............................

c     null the species number densities and ion ionization fractions

      DO 51 k=1,Imax
       nk(k)       = 0.0d0
       alphksol(k) = 0.0d0
       fabund(k)   = 0.0d0
       xfrac(k)    = 0.0d0
       DO 52 j=1,Imax+1
        nkj(k,j)  = 0.0d0
        fion(k,j) = 0.0d0
 52    CONTINUE
 51   CONTINUE

c     ...........
c     TIME SCALES [for target ions kdo(m),jdo(m)]
c     ...........

      DO 61 m=1,Imax
       tau_ph(m)   = 0.0d0     ! photoionization timescale
       tau_rec(m)  = 0.0d0     ! recombination timescale
       tau_coll(m) = 0.0d0     ! collisional iomnization timescale
 61   CONTINUE

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE sappend(instring,appstring,outstring)

c     append appstring onto instring and output the result as outstring

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             i,k,lend
      character*(*)       instring,appstring,outstring

      lend = len(instring)
      k    = 0
      DO 09 i=1,lend
        k = i
        IF (instring(i:i).eq.' ') GOTO 10
 09   CONTINUE

 10   outstring = instring(1:k-1)//appstring

      RETURN
      END


c
c
c.........................................................................
c

      SUBROUTINE notify(instring,n)

c     communicate to screen and log file 

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             n,lens
      character*(*)       instring

      lens = len(instring)

      IF (n.eq.-1) then
       WRITE(6,*) instring(1:lens)
       WRITE(4,*) instring(1:lens)
      ELSE
       WRITE(6,*) instring(1:lens),n
       WRITE(4,*) instring(1:lens),n
      ENDIF

      RETURN
      END


c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     the Numerical Recipes for interpolation and integration
c
c     interpolations routines:
c     SUBROUTINE SPLINE
c     SUBROUTINE SPLINT
c
c     integration routines:
c     SUBROUTINE QROMB
c     SUBROUTINE TRAPZD
c     SUBROUTINE POLINT
c
c     root solve routines
c     FUNCTION ZBRENT
c
c     the only mods are 
c     (1) implicit double precision declarations
c     (2) error messages and stops
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

c
c.........................................................................
c

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)

c     computes 1D 2nd derivatives of array Y with respect to array X

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=10000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)

      IF (YP1.GT..99D30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF

      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        IF (SIG.le.0.0) GOTO 99
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE

      IF (YPN.GT..99D30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)

      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE

      RETURN

 99   WRITE(6,*) 'ERROR(spline): X data ill-posed index I=',I
      WRITE(6,*) 'X data must be in ascending order'
      STOP

      END

c
c.........................................................................
c

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)

c     performs 1D interpolation of YA given X; returns Y

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(N),YA(N),Y2A(N)

      KLO=1
      KHI=N

1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF (XA(K).GT.X) THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

      H=XA(KHI)-XA(KLO)

      IF (H.EQ.0.) THEN
       WRITE(6,*) 'ERROR(splint): X cannot be bracketed due to'
       WRITE(6,*) 'bad XA input at indices KLO,KHI= ',KLO,KHI
       STOP
      END IF

      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.

      RETURN
      END



c
c.........................................................................
c

      SUBROUTINE QROMB(FUNC,A,B,SS)

c     performs Romberg's interative refinement integration of function
c     FUNC from A to B; returns the refined integral as S; we use the
c     routine TRAPZD (see below)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(EPS=1.E-7,JMAX=50,JMAXP=JMAX+1,K=5,KM=4)
      DIMENSION S(JMAXP),H(JMAXP)

      H(1)=1.

      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
	  L=J-KM
          CALL POLINT(H(L),S(L),K,0.,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25*H(J)
11    CONTINUE

      WRITE(6,*) 'ERROR(qromb): too many steps trying to'
      WRITE(6,*) 'converge the integral to a tolerance of'
      WRITE(6,*) eps
      WRITE(6,*) 'Assuming the function is not ill-posed,'
      WRITE(6,*) 'this can be remedied by decreasing EPS'
      WRITE(6,*) 'and/or increasing JMAX in routine QROMB'
      STOP 

      END


c
c.........................................................................
c

      SUBROUTINE TRAPZD(FUNC,A,B,S,N)

c     performs trapezoid rule integration of function FUNC from A to B;
c     returns the integral as S

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE IT

      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF

      RETURN
      END


c
c.........................................................................
c

      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)

c     performs polynomial interpolation to estimate the error in the
c     integral computed in routine QROMB

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=10) 
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)

      NS=1
      DIF=ABS(X-XA(1))

      DO 11 I=1,N 
        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I
          DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
11    CONTINUE

      Y=YA(NS)
      NS=NS-1

      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          IF (DEN.EQ.0.) GOTO 99
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
12      CONTINUE
        IF (2*NS.LT.N-M) THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
13    CONTINUE

      RETURN

 99   WRITE(6,*) 'ERROR(polint): division by zero error.'
      STOP 

      END

c
c.........................................................................
c

      FUNCTION ZBRENT(FUNC,X1,X2,TOL)

c     root solves the function FUNC between the points X1 and X2 to a
c     tolerance of TOL

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ITMAX=500,EPS=3.E-16)

      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)

      IF (FB*FA.GT.0.) GOTO 99

      FC=FB

      DO 11 ITER=1,ITMAX

        IF (FB*FC.GT.0.0) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF

        IF (ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=C
          C=A
          FA=FB
          FB=FC
          FC=FA
        ENDIF

        TOL1=2.0d0*EPS*ABS(B)+0.5d0*TOL
        XM=0.5d0*(C-B)

        IF (ABS(XM).LE.TOL1 .OR. FB.EQ.0.0) THEN  ! good root
          ZBRENT=B
          RETURN
        ENDIF

        IF (ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.0d0*XM*S
            Q=1.0d0-S
          ELSE
            Q=FA/FC
            R=FB/FC
            P=S*(2.0d0*XM*Q*(Q-R)-(B-A)*(R-1.0d0))
            Q=(Q-1.0d0)*(R-1.0d0)*(S-1.0d0)
          ENDIF
          IF (P.GT.0.0) Q=-Q
          P=ABS(P)
          IF (2.0d0*P .LT. MIN(3.0d0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF

        A=B
        FA=FB

        IF (ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF

        FB=FUNC(B)

11    CONTINUE

c     iterations exceeded, bad root

      XM     = 0.5d0*(C-B)
      TOL1   = 2.0d0*EPS*ABS(B)+0.5d0*TOL
      ZBRENT = 0.5d0*(B+C)

      WRITE(6,*) 'WARNING(zbrent): exceeded maximum iterations-'
      WRITE(6,*) 'returning the average of the last bracketing'
      WRITE(6,*) 'values B,C     = ',B,C
      WRITE(6,*) 'where (C-B)/2  = ',XM
      WRITE(6,*) 'and    TOL1    = ',TOL1
      WRITE(6,*) 'and    TOLER   = ',TOL
      WRITE(6,*) 'and returned root = ',zbrent
      WRITE(6,*) 'Root should be viewed with skepticism.'
      RETURN

 99   WRITE(6,*) 'ERROR(zbrent): Root must be bracketed for ZBRENT.'
      WRITE(6,*) 'The bracket value (X1 or X2) that is closest to'

      IF (abs(func(X1)).lt.abs(func(x2))) then
       WRITE(6,*) 'the root is X1,F(X1) = ',x1,func(x1)
       WRITE(6,*) 'where X2,F(X2)       = ',x2,func(x2)
      ELSE
       WRITE(6,*) 'the root is X2,F(X2) = ',x2,func(x2)
       WRITE(6,*) 'where X1,F(X1)       = ',x1,func(x1)
      END IF
      STOP

      END

c..............................................................................
c

      SUBROUTINE          today(adat,atim)   

c
c     this routine gets the date and time out of a machine.  the output
c     format is

c     adat dd mmmyy
c     atim hh:mm:ss

c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      character*3        amon(12) 
      character*12       atim,adat
      character*10       adigit   
      parameter          (adigit='0123456789')
      integer*4          idat(3),itim(3)  

      data amon/ 'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,    
     +           'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec' /

      adat=' '  
      atim=' '  

c     grab the time and "write" into variable atim

      CALL itime(itim)  
      write (atim,113) itim  

c     grab the date and write into variable adat

      CALL idate(idat)  
      write (adat,114) idat(1),amon(idat(2)),idat(3)

 113  format(i2,':',i2,':',i2)  
 114  format(i2,'-',a3,' ',i4)  

      return

      end   

