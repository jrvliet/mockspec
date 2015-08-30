c
c.........................................................................
c

      SUBROUTINE Rkjsgrid(ilev,kdo,jdo,toler,oflag)

c     for cell level ilev, this routine sets up the nH and T axes for
c     the grid of ionization fractions for the target ions

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      include           'getcube.h'
      logical           oflag
      integer           ilev,in,iT,itau,j,k,s,kk,nek
      integer           kdo(Imax),jdo(Imax)
      double precision  toler,T,Zmet
      character*12       adate,atime
      include           'rates.com'
      include           'com-photo.com'
      include           'com-phxsecs.com'
      include           'com-auger.com'
      include           'com-Rkjs.com'
      include           'getcube.com'


c     initialize the boolean that flags if there are optically thick
c     cells for this cell level; if OFLAG remains false then no
c     self-shielding grid is required for the current cell level and we
c     return to the calling routine

      oflag = .false.

c     populate the photoionization and Auger rates for the unattenuated
c     SED

      CALL pop_photoRs

c     set up the grid for this level; routine mknHTgrid checks all the
c     cells one by one to see if any of them are optically thick; if so,
c     then the grid is built; OFLAG is also determined in mknHTgrid

      CALL zerofgrid
      CALL mknHTgrid(ilev,Zmet,toler,oflag)

c     are their optically thick cells on this level?  If not, then
c     return

      IF (.not.oflag) RETURN

c     ok, so we have optically thick cells, build the grid 
c     [this is time consuming!]

      CALL notify('  constructing grid...',-1)
      CALL mkfgrid(ilev,kdo,jdo,Zmet,toler)
      CALL notify('  grid complete, solving cells...',-1)

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE mknHTgrid(ilev,Zmet,toler,oflag)

c     this sets up the nH and T axes for the ionization fraction grid
c     for cells of level ILEV

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      include           'getcube.h'
      logical           opticallythick,oflag
      integer           ilev,i,in,iT
      double precision  Zmet,toler
      double precision  dnH,dT,nHlo,nHhi,Tlo,Thi,T,pow
      include           'rates.com'
      include           'getcube.com'
      include           'com-Rkjs.com'


c     obtain the geometric mean metallicity for the grid of this level
c     TEMPORARY APPROACH (I have not tested how sensitive the helium
c     edges are to metal mass fraction; however, because the total mass
c     fraction is conserved, as Zmet decreases, the abundance of helium
c     goes up, so lower Zmet grids will have higher helium edge opacity
c     and visa versa)

      Zmet = sqrt(Zmin(ilev)*Zmax(ilev))

c     STEP 1: search phase space with 0.1 log10 intervals for the
c     mininum and maximum nH and T that are optically thick for this
c     level

      dnH  = 0.1
      dT   = 0.1      
      NRn  = int(log10(nHmax(ilev)/nHmin(ilev))/dnH) + 1
      NRT  = int(log10(Tmax(ilev)/Tmin(ilev))/dT) + 1
      nHlo = nHmax(ilev)
      nHhi = nHmin(ilev)
      Tlo  = Tmax(ilev)
      Thi  = Tmin(ilev)

      DO in=1,NRn          
       pow  = float(in-1)*dnH
       nH   = nHmin(ilev)*(10.0d0**pow)         
       CALL initabund(Zmet)
       DO iT=1,NRT           
        pow  = float(iT-1)*dT
        T    = Tmin(ilev)*(10.0d0**pow)         
        CALL alphas_T(T)     
        CALL ionbalance(toler) 
        IF (opticallythick(ilev)) then ! check if optically thick
         nHlo  = min(nHlo,nH)
         nHhi  = max(nHhi,nH)
         Tlo   = min(Tlo,T)
         Thi   = max(Thi,T)
         oflag = .true.
        END IF
       ENDDO
      ENDDO

c     if none of the cells are optically thick, then communicate and
c     return

      IF (.not.oflag) then
       CALL notify(' all cells optically thin',-1)
       RETURN
      ENDIF

c     STEP 2: create the final grid

c     set the nH and T boundary values for the density and temperature
c     axes of the grid, these are now stored as log10; we add an
c     additional increment in the ranges to be sure to bracket all cells

      IF (nHlo.gt.nHmin(ilev)) then          ! minimum nH
       nHmin(ilev) = log10(nHlo)-dnH
      ELSE
       nHmin(ilev) = log10(nHmin(ilev))
      ENDIF

      IF (nHhi.lt.nHmax(ilev)) then          ! maximum nH
       nHmax(ilev) = log10(nHhi)+dnH
      ELSE
       nHmax(ilev) = log10(nHmax(ilev))
      ENDIF

      IF (Tlo.gt.Tmin(ilev)) then            ! minimum T
       Tmin(ilev)  = log10(Tlo)-dT
      ELSE
       Tmin(ilev)  = log10(Tmin(ilev))
      ENDIF

      IF (Thi.lt.Tmax(ilev)) then            ! maximum T 
       Tmax(ilev)  = log10(Thi)+dT
      ELSE
       Tmax(ilev)  = log10(Tmax(ilev))
      ENDIF

c     set up the nH axis; the number of elements and values

      NRn  = int((nHmax(ilev)-nHmin(ilev))/dnH) + 1
      IF (Nrn.gt.Nrnmax) then
       WRITE(6,*) 'ERROR(mkgrid): Number of nH gird points ',NRn
       WRITE(6,*) 'exceeds maximum allotment NRnmax = ',NRnmax
       WRITE(6,*) 'To fix, increase NRnmax and recompile.'
       STOP
      ENDIF
      DO i=1,NRn
       nR(i) = nHmin(ilev) + float(i-1)*dnH
      ENDDO

c     set up the T axis; the number of elements and values

      NRT  = int((Tmax(ilev)-Tmin(ilev))/dT) + 1
      IF (Nrn.gt.Nrnmax) then
       WRITE(6,*) 'ERROR(mkgrid): Number of T gird points ',NRT
       WRITE(6,*) 'exceeds maximum allotment NRTmax = ',NRTmax
       WRITE(6,*) 'To fix, increase NRnmax and recompile.'
       STOP
      ENDIF
      DO i=1,NRT
       TR(i) = Tmin(ilev) + float(i-1)*dT
      ENDDO

c     the opacity grid [this is inefficient to place here, I should move
c     it later so that it is called only once in the whole program]

      tau(1) = 1.0d0
      tau(2) = 1.0d0
      tau(3) = 1.0d0
      tau(4) = 1.0d0
      tau(5) = 1.0d0
      DO i=6,NRtaumax
       IF (tau(i-1).lt.100.0) then
         tau(i) = 1.30d0*tau(i-1)
       ELSE
         tau(i) = tau(i-1)
       END IF
      ENDDO

c     communicate the grid structure

      CALL notify('  sheilding grid structure is',-1)
      WRITE(6,600) NRn,nHmin(ilev),nHmax(ilev),dnH
      WRITE(6,601) NRT,Tmin(ilev),Tmax(ilev),dT
      WRITE(6,602) 1,log10(Zmin(ilev)),log10(Zmax(ilev)),log10(Zmet)
      WRITE(4,600) NRn,nHmin(ilev),nHmax(ilev),dnH
      WRITE(4,601) NRT,Tmin(ilev),Tmax(ilev),dT
      WRITE(4,602) 1,log10(Zmin(ilev)),log10(Zmax(ilev)),log10(Zmet)

      RETURN

 600  FORMAT(4x,'log(nH) axis: N =',i3,' range = (',f6.3,',',f6.3,')',
     &          '  step = ',f6.3)
 601  FORMAT(4x,'log(T)  axis: N =',i3,' range = (',f6.3,',',f6.3,')',
     &          '  step = ',f6.3)
 602  FORMAT(4x,'log(Zmet)   : N =',i3,' range = (',f6.3,',',f6.3,')',
     &          '  Zmet = ',f6.3)

      END

c
c.........................................................................
c

      SUBROUTINE save_photoRs

c     this routine populates the partial (i.e., shell) photoionization
c     rates for ion k,j for optical depth zero (face of the cell); this
c     routine is called once after the SED is built and before the data
c     cube is read

c     these rates depend only on the shape of the ionizing SED

c     if using stellar SEDs too, this will require modifying

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      include           'getcube.h'
      integer            j,k,s,kk,nek
      double precision  PH_Rate
      include           'rates.com'
      include           'com-Rkjs.com'

c     we need the total number of populated shells, so include the ntot

      integer           ntot
      COMMON/ntot/      ntot(30)


c     fill the grid for the unattenuated SED; the function PH_rate
c     computes and stores the partial photoionization rates, R_phs(s),
c     as a global variable at each k,j; we store both the total
c     photoionization rate, R_ph, for ion k,j, and the partial rates,
c     Rphkjs for each shell s of ion k,j

c     outer loop is over species

      DO kk=1,Nspecies

       k = kidx(kk)

c     loop over ionization stage; evaluate the rates (function PH_Rate)

        DO j=1,k

         R_ph(k,j) = PH_Rate(k,j)

c     compute the number of electrons for the iselectronic sequence of
c     ion k,j and loop over the shells and store the partial
c     photoionization rates, which were compute in function PH_Rate and
c     stored globally

         nek = k - j + 1
         DO s=1,ntot(nek)
          Rphkjs0(k,j,s) = R_phs(s)
        ENDDO     

       ENDDO      ! next ionization stage
      ENDDO       ! next species

      RETURN
      END
c
c.........................................................................
c

      SUBROUTINE pop_photoRs

c     this routine populates the photoionization and Auger rates for ion
c     k,j for optical depth zero (face of the cell); basically we unfold
c     the data stored in routine save_photoRs

c     once we populate the rates, we can solve optically thin cells or
c     the facial layer of optically thick cells

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      include           'getcube.h'
      integer            j,k,s,kk,nek
      double precision   AUG_Rate
      include           'rates.com'
      include           'com-Rkjs.com'
      include           'com-auger.com'
      include           'com-modes.com'

c     we need the total number of populated shells, so include the ntot

      integer           ntot
      COMMON/ntot/      ntot(30)

c     loop over species, ionization stage, and shell and populate the
c     partial and photoionization rates; then use the partial rates to
c     compute the Auger rates, Q and R_Agrout

c     outer loop is over species; obtain the index of the species; then
c     loop over ionization stage; the photoionization rate for species
c     j,k for the summation over the partial (shell) rates

      DO kk=1,Nspecies

       k = kidx(kk)

        DO j=1,k

c     loop over the shells, store the partial rates, R_phs, and compute
c     the total rate, R_ph

         nek = k - j + 1
         DO s=1,ntot(nek)          
          R_phs(s)  = Rphkjs0(k,j,s) 
          R_ph(k,j) = R_ph(k,j) + W_Auger(k,j,s,1)*R_phs(s)
         ENDDO     

c     if doing Auger ionization, compute the rates

        IF (doAUG) then
          CALL getQs(k,j)
          R_Agrout(k,j) = AUG_Rate(k,j)
        ENDIF

       ENDDO       ! next ionization stage

      ENDDO        ! next species

      RETURN
      END


c
c.........................................................................
c
      SUBROUTINE  zerofgrid

c     this routine is called each time we start a new cell level; it's
c     puprose is to initialize/null the grid used to treat self
c     shielding

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none     
      include           'rates.h'
      integer           in,iT,m
      include           'com-Rkjs.com'

      DO in=1,NRnmax 
       nR(in) = 0.0d0                         
       DO iT=1,NRTmax                         
        TR(iT) = 0.0d0                         
        DO m=1,Imax
         fgrid(m,in,iT) = 0.0d0
        ENDDO                      
       ENDDO                               
      ENDDO                                

      RETURN
      END
