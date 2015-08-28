c
c	Subroutines from ART
c--------------------------------------------------
CEVERINO07052005
c	I/O SUBROUTINES
C----------------------------
C
c     --------------------------------------------
      subroutine Read_Gas_Binary ( dt , FileName )
c     --------------------------------------------
c
c     purpose: reads files with control information and gas data
c              and hydrodynamic variables for cells
CEVERINO07042005: Disable previous time moment for the energy cons. control
c

      include 'a_tree.h'
      include 'a_control.h'
      real*8 dt
      character*256 FileName
      character*256 jname
      real*8 b2a

      nfname = index ( FileName , ' ' ) - 1
      open  ( 19 , file = FileName(1:nfname) , form = 'unformatted' )
      write(*,*) 'reading tree data from file ',FileName(1:nfname)
      read ( 19 ) jname 
      if ( jname .ne. jobname1 ) then 
        write(*,*) '* warning : Read_Gas_Binary : jname .ne. jobname1'
        write(*,'(2(A256))') jname, jobname1
c        stop
      endif
      read ( 19 ) istep , t , dt, aexpn, ainit  
      write(*,20) 
     &      istep, t, dt, aexpn, ainit
 20   format ('istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)
      do Level = MinLevel , MaxLevel 
        aexp(Level) = aexpn
      enddo
c     if having trouble on linux, try inserting a dummy first in the read
c     list below
      read ( 19 )  boxh, Om0, Oml0, Omb0, hubble
      write(*,21) 
     &      boxh, Om0, Oml0, Omb0, hubble      
 21   format('Lbox[/h Mpc] =',f9.4,1x,' Om0 = ',f5.3,1x,
     &        ' Oml0= ',f5.3,1x,' Omb0= ',f7.5,1x,' h = ',f5.3)
      read ( 19 ) nextras
      if ( nextras .gt. nextra ) then 
        write(*,*) '* error: Read_Gas_Binary : nextras > nextra :',
     &             nextras, nextra
        write(*,*) '* (check nextra parameter in a_control.h) '
        stop
      endif
      read ( 19 ) (extra(i), i=1,nextras)
      read ( 19 ) (lextra(i), i=1,nextras)

      read ( 19 ) MinLev , MaxLevelNow
      if ( MaxLevelNow .gt. MaxLevel ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : MaxLevelNow > MaxLevel :',
     &        MaxLevelNow, MaxLevel 
         stop
      endif
      read ( 19 ) ( tl(i)     , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtl(i)    , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( tlold(i)  , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( dtlold(i) , i = MinLevel , MaxLevelNow )
      read ( 19 ) ( iSO(i)    , i = MinLevel , MaxLevelNow )
      read ( 19 ) ncell
      if ( ncell .ne. ncell0 ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : ncell not equal ncell0 :',
     &        ncell, ncell0
         stop        
      endif
      read ( 19 ) (iOctCh(i), i=1,ncell)
      read ( 19 ) ((hvar(ivar,i),ivar=1,nhvar), i=1,ncell)
      read ( 19 ) ((var(ivar,i),ivar=2,3), i=1,ncell)
      ntot = 0 
      if ( MaxLevelNow .gt. MinLevel ) then 
        read ( 19 ) iOctFree , nOct
        write(*,*) 'iOctFree =',iOctFree,' nOct =', nOct
        if ( ncell+nchild*nOct .gt. mcell ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : ncell+nchild*nOct > mcell :',
     &        ncell+nchild*nOct, mcell 
         write(*,*)  
     &    '* (size of the input exceeds max. of cells (mcell))'
         write(*,*) '* increase mcell in a_setup.h and try again...'
         stop
        endif 

        do Lev = MinLevel+1 , MaxLevelNow
          read ( 19 ) Level , iNOLL(Level) , iHOLL(Level)
          iOct   = iHOLL(Level)
          nLevel = iNOLL(Level)
          nLevCells = nLevel*nchild
          if ( nLevCells .gt. nclmax ) then 
            write(*,*) 'error : L =',Lev,' nLevCell =',nLevCell,
     &        ' > nclmax =',nclmax,' set in a_setup.h'
            write(*,*) '=> increase nclmax and rerun.'
            close ( 19 )
            stop
          endif
          ntot = ntot + nLevel 
          write(*,*) 'reading tree (oct) data for level ', Lev
          do ic1 = 1 , nLevel
            read(19) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
     &                iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct), 
     &                iOctLL2(iOct)
            iOct = iOctLL1(iOct)
          enddo
          write(*,*) 'reading cell data: ncells =', nLevel*nchild
          do ic1 = 1 , nLevel*nchild
            read ( 19 ) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar), 
     &                 (var(i,idc), i=2,3)
          enddo
        enddo
      endif
      close (19)
      write(*,*) 'done reading',ncell0+nchild*ntot,' cells...'

c     set previous time moment for the energy cons. control
c      ap0    = b2a ( (1.d0*(t-5.d-1*dt)) )      
c
      return
      end
c
c     ----------------------------------------------------------
      subroutine Read_Particles_Binary ( FileName1 , FileName2 )
c     ----------------------------------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
c
      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
      character*256 FileName1 , FileName2
c 
      ngrid = ng

      nfn1 = index ( FileName1 , ' ' ) - 1
      nfn2 = index ( FileName2 , ' ' ) - 1
c
      write (*,*)  'Open file...', FileName1(1:nfn1)
      open ( 3 , file = FileName1(1:nfn1) , form = 'unformatted' )
      write (*,*)   FileName1(1:nfn1), '  is open.'
c.... read control information and check whether it has proper structure

      read      (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,hubble,Wp5,
     &                  Ocurv,Omb0,extras
      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRIDC,NRECL,Om0,Oml0,Omb0,hubble
 100  format (1X,'Header=>',A45,/
     &           1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &           1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &           1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     &           1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,
     &               ' Omegab_0=',F7.3,' Hubble=',f7.3)

      write(*,*)'extras: specie   wspecie    lspecie '
      do ispec=1,nspecies
         write (*,*)  ispec, wspecies(ispec) , lspecies(ispec)    
      enddo
c...  Working with new ICs (Daisuke added a line below) 
      if( nspecies .eq. 0 ) nspecies = 1
      if ( nspecies .gt. nspec ) then
        write (*,*)
     &      ' nspec in PARAMETER < nspecies in TAPE-FILE:'
         write (*,*) ' nspecies=',nspecies,' > ',nspec
         STOP
      endif
      If( npmax .lt. lspecies(nspecies) ) then
         write (*,*) ' Wrong number of particles: '
         write (*,*) ' should be =',lspecies(nspecies),' (lspecies)'
         write (*,*) ' but is set to ',npmax,' in a_setup.h...'
         do ispec = 1 , nspecies
           write(*,*) ispec, lspecies(ispec)
         enddo
         STOP
      Endif 

      nbyte  = nrecl * 4
      nacces = nbyte
#ifdef NBYTEWORD
      nacces = nrecl
#endif

      open ( 1 , file = FileName2(1:nfn2), access = 'direct',
     &	         status = 'unknown', recl = nacces      )
 
      rewind 3

      N_particles = lspecies(nspecies)
      np = N_particles 
      write(*,*) ' hello :',N_particles, np
      do i = 1 , nspecies
        write(*,*) 'ispecie =',i,lspecies(i)
      enddo
      write(*,*) 'Read_Particles_Binary: reading np=',np,
     &                ' particles'


      if ( np .eq. 0 .or. np .gt. npmax ) then 
        write(*,*) '1 : error in Read_Particles_Binary: bad np =',
     &             np
        write(*,*) 'stopping...'
        stop
      endif

      Npages      = (N_particles -1)/npage + 1
      N_in_last   = N_particles - npage*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',nspecies
      write (*,*) ' N_in_last=',N_in_last

      do irow = 1 , Npages         ! loop over particle pages
        In_page = npage
        if ( irow .eq. Npages ) In_page = N_in_last
c         write (*,*)' Read page=',IROW,' file=',ifile,' N=',In_page
        iL = npage * (irow-1)
        CALL GetRow(irow,1) ! read in a page of particles
        do in = 1 , In_page          ! Loop over particles
          ip = in + iL                     ! current particle number
          x(ip) = xpar(in)
          y(ip) = ypar(in)
          z(ip) = zpar(in)
          vx(ip) = vxx(in)
          vy(ip) = vyy(in)
          vz(ip) = vzz(in)
        enddo
      enddo

      close (1)
      close (3)

      do i = 1 , nspecies
        wpar(i) = wspecies(i)
        lsp(i)  = lspecies(i)
      enddo
      do i = 1 , nspecies
        if(i .eq. 1 ) then
          nsp(i,1) = 1
        else
           nsp(i,1) = lsp(i-1) + 1
        endif 
        nsp(i,2) = lsp(i)
      enddo 

C$OMP PARALLEL DO 
      do ic1 = 1 , np 
        pt(ic1) = t   ! set t to the current time moment
      enddo
c
c.... compute each particle's weight
c
      n_spec = nspecies
#ifdef STARFORM
      n_spec = nspecies-1
#endif
      do ic0 = 1 , n_spec 
        do ic1 = nsp(ic0,1) , nsp(ic0,2)
           pw(ic1) = wpar(ic0)
        enddo
      enddo 

c$$$#ifdef STARFORM
c$$$C      if ( start ) then         
c$$$        if ( nspecies+1 .gt. nspec ) then
c$$$          write (*,*)
c$$$     &        ' nspec in PARAMETER < nspecies+1 for starformation:'
c$$$          write (*,*) ' nspecies+1=',nspecies+1,' > nspec =',nspec
c$$$          STOP
c$$$        endif
c$$$        nspecies = nspecies + 1
c$$$        nsp(nspecies,1) = lsp(nspecies-1) + 1
c$$$        nsp(nspecies,2) = lsp(nspecies-1)
c$$$        lsp(nspecies) = nsp(nspecies,2)
c$$$        lspecies(nspecies) = lsp(nspecies)
c$$$
c$$$
c$$$        write(*,*) 'in reading particles:'
c$$$        write(*,*) '    lsp      wpar         nsp1      nsp2 '
c$$$        do i = 1 , nspecies
c$$$           write(*,30) lsp(i), wpar(i), nsp(i,1), nsp(i,2)
c$$$        enddo
c$$$        write(*,*) ' '
c$$$C      endif
c$$$#endif
c$$$      
c$$$ 30   format (i10,2x,g11.5,2x,2(i8,2x))

      return
      end
c     
c     ------------------------------------------
      subroutine Read_Stellar_Data ( FileName )
c     ------------------------------------------
c
c     purpose: when the code starts or restarts
c              this routine reads all data on stellar particles 
c 
      include 'a_tree.h'
      include 'a_control.h'
      character*256 FileName 
      real*8 tdum, adum 
      real*8 ws_old, ws_oldi
      common / s_old1 / ws_old, ws_oldi
c
      nfn = index ( FileName , ' ' ) - 1
      open ( 60 ,file = FileName(1:nfn),
     &       form = 'unformatted',status = 'old' )
      read(60) tdum, adum
      read(60) nstars

      if ( nstars .ne. nsp(nspecies,2) - nsp(nspecies,1) + 1 ) then 
        call Open_ASCII_File ( 13 , errorfile , sequent ) 
        write(13,*) '1  : In Read_Stellar_Data: something is wrong:'
        write(13,*) 'nstars =',nstars,' is iconsistent with'
        write(13,*) 'nsp(nspecies,2)-nsp(nspecies,1)+1 =',
     &               nsp(nspecies,2) - nsp(nspecies,1) + 1
        write(13,*) 'stopping...'
        close(13)              
        stop
      endif
      if ( nstars .eq. 0 ) return
      if ( nstars .lt. 0 .or. nstars .gt. nstarmax ) then 
        call Open_ASCII_File ( 13 , errorfile , sequent ) 
        write(13,*) '1  : In Read_Stellar_Data: bad nstars =',nstars
        write(13,*) 'must be 0 < nstars < nstarmax =',nstarmax        
        write(13,*) 'stopping...'
        close(13)              
        stop
      endif
      read (60) ws_old, ws_oldi
      read (60) (pw(ic1),ic1=nsp(nspecies,1),nsp(nspecies,2))   ! weights
      read (60) (pw0(ic1),ic1=1,nstars)      ! initial masses
      read (60) (tbirth(ic1),ic1=1,nstars)   ! birth times
#ifdef ENRICH
      read (60) (zstII(ic1),ic1=1,nstars)    ! metallicity 
#endif
#ifdef ENRICH_SNIa
      read (60) (zstIa(ic1),ic1=1,nstars)    ! metallicity 
#endif
      close( 60 )
c
      return
      end
c
c     ------------------------------------------
      subroutine Write_Gas_Binary ( dt , fname )
c     ------------------------------------------
c
c     purpose: writes files with control information and gas data
c              writes hydrodynamic variables for cells
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      character*256 fname
      real*8 b2a
      real*8 dt
c
      nfn = index ( fname , ' ' ) - 1
      call Get_MaxLevelNow ()
      open  ( 22 , file = fname(1:nfn) , form = 'unformatted' )
      write ( 22 ) jobname1
#ifdef COSMOLOGY
      adum = b2a ( (1.d0*t) )
#endif
#ifndef COSMOLOGY
      adum = 1.0 
#endif
      write ( 22 ) istep , t , dt , adum , ainit
      write ( 22 ) boxh, Om0, Oml0, Omb0, hubble
      write ( 22 ) nextra
      write ( 22 ) (extra(i), i=1,nextra)
      write ( 22 ) (lextra(i), i=1,nextra)
      write ( 22 ) MinLevel , MaxLevelNow
      write ( 22 ) ( tl(i)     , i = MinLevel , MaxLevelNow )
      write ( 22 ) ( dtl(i)    , i = MinLevel , MaxLevelNow )
      write ( 22 ) ( tlold(i)  , i = MinLevel , MaxLevelNow )
      write ( 22 ) ( dtlold(i) , i = MinLevel , MaxLevelNow )
      write ( 22 ) ( iSO(i)    , i = MinLevel , MaxLevelNow )
      write ( 22 ) ncell0
      write ( 22 ) (iOctCh(i), i=1,ncell0)
      write ( 22 ) ((hvar(ivar,i),ivar=1,nhvar), i=1,ncell0)
      write ( 22 ) ((var(ivar,i),ivar=2,3), i=1,ncell0)
      if ( MaxLevelNow .gt. MinLevel ) then 
        write ( 22 ) iOctFree , nOct
        do Level = MinLevel+1 , MaxLevelNow
          nLevel = iNOLL(Level)      
          write ( 22 ) Level , iNOLL(Level) , iHOLL(Level)
          iOct = iHOLL(Level)
          do ic1 = 1 , nLevel
            write(22) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
     &                iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct), 
     &                iOctLL2(iOct)
            iOct = iOctLL1(iOct)
          enddo
          call Select_Cells ( Level , nLevel ) 
          do ic1 = 1 , nLevel 
            icell = iSelect(ic1)
            do ic2 = 0 , 7
              idc = icell + ic2
              write(22) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar), 
     &                 (var(i,idc), i=2,3)
            enddo
          enddo
        enddo
      endif
      close (22)
      write(*,*) 'done writing...'
c
      return
      end
c
c     --------------------------
      real*8 function b2a ( bt )
c     --------------------------
c
      real*8 bt
      real*8 f_b2a, fp(1)
      real zbrent 
      external zbrent, f_b2a
c
      include 'a_numbers.h'
      include 'a_setup.h'
      include 'a_control.h'
c
      IF ( (Om0 .eq. one) .and. (Oml0 .eq. zero) ) THEN 
        b2a = (1.d0 / (1.d0 - bt))**2
      ELSE
        if ( (bt .lt. -190.d0) .or. (bt .gt. 1.d-1) ) then 
          write(*,*) ' b2a : something is wrong : bt =',bt
          stop
        endif
        if ( (Om0 .eq. zero) .and. ( Oml0 .eq. zero ) ) then 
          write(*,*) 'b2a : something is wrong : Om0 = Oml0 = 0.0'
          stop
        endif
        fp(1) = bt
        b2a = zbrent ( f_b2a , fp , 1 , 1.d-4 , 1.1d0 , 1.d-9 )
      ENDIF
c
      return
      end
c
c     --------------------------------------
      real*8 function f_b2a ( at , fp , np )
c     --------------------------------------
c
c     at - expansion factor
c     fp(1) = b
c
      integer np
      real*8 at , fp(np)
      real*8 a2b
      external a2b
c
      f_b2a = a2b ( at ) - fp(1)
c      
      return
      end
c
c     --------------------------
      real*8 function a2b ( at )
c     --------------------------
c
c     translates expansion factor at into hydro time variable b
c
      real*8 at
      real*8 fp(2), atst
      real*8 INTEGRATE , f_a2b
      external INTEGRATE , f_a2b

      include 'a_numbers.h'
      include 'a_setup.h'
      include 'a_control.h'
c      
c
      if ( (at .lt. zero) ) then 
        write(*,*) 'a2b : something is wrong : at = ', at
        stop
      endif
      IF ( (Om0 .eq. one) .and. (Oml0 .eq. zero) ) THEN 
        a2b = (1.d0 - 1.d0/sqrt(at))
      ELSE
        if ( (Om0 .eq. zero) .and. ( Oml0 .eq. zero ) ) then 
          write(*,*) 'a2b : something is wrong : Om0 = Oml0 = 0.0'
          stop
        endif
        fp(1) = Om0
        fp(2) = Oml0
        atst  = 1.d-1 * (1.d0 - at)
        a2b = INTEGRATE( f_a2b , fp , 2 , 1.d0 , at , atst , 1.d-9 )
      ENDIF       
c
      return
      end
c
c     -------------------------------------
      real*8 function f_a2b ( x , fp , np )
c     -------------------------------------
c
c     input : x     - is expansion factor variable
c             fp(1) - Om0  = present-day matter density 
c             fp(2) - Oml0 = present-day vacuum contribution
c
      integer np 
      real*8 x , fp(np), d1
c
      d1 = x**3
      f_a2b = .5d0 * sqrt(fp(1)) / d1 
     &             / sqrt(fp(1)/d1 + fp(2) + (1.d0-fp(1)-fp(2))/x**2)
c
      return
      end
c
C-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE units()
c-------------------------------------------------------------------- 
c     set units for physical variables
c
#     include "a_def.h"
      include 'a_numbers.h'
      include 'a_setup.h'
      include 'a_control.h'
      real wmu53

      wmu53 = wmu**(5./3.)

      if ( boxh .le. 0. ) then 
        write(*,*) '* SetUnits : something is wrong: boxh <= 0 :',boxh
        stop
      endif

c
c.... comoving size of the zeroth-level cell in /h Mpc
c
      r0 = boxh / ng 

      if ( Om0 .le. 0. ) then 
        write(*,*) '* SetUnits : something is wrong: Om0 <= 0 :',Om0
        stop
      endif
      if ( hubble .le. 0. ) then 
        write(*,*)'* SetUnits : something is wrong: hubble <= 0:',hubble
        stop
      endif
      if ( hubble .gt. 1. ) then 
        write(*,*)'* warning : SetUnits : suspicious hubble =',hubble
      endif

      dummy = 100.0 * hubble * sqrt(Om0) 

      
c.... time unit in yrs
c 
      t0 = 2.0 / dummy * 3.0856e19 / 3.15e7
c
c.... v0 - velocity units in km/s
c
      v0 = 50. * r0 * sqrt(Om0)
c
c.... rho0 = 3H0^2 * Om0 / (8*pi*G) - unit of density in Msun/Mpc^3
c
      rho0 = 2.776e11 * hubble**2 * Om0
c
c.... P0 = rho0 * v0**2 - unit of pressure in g/cm/s^2
c
      P0 = 4.697e-16 * Om0**2 * r0**2 * hubble**2
c
c.... T_0 = unit of temperature in K and in keV)
c
      T_0 = 3.03e5 * r0**2 * wmu * Om0 ! [K]
c      T_0 = 2.61155 * r0**2 * wmu * Om0 ! [keV]
c
c.... S_0 = unit of entropy in keV * cm^2
c
     
      S_0 = 52.077 * wmu53 * hubble**(-4./3.)*Om0**(1./3.)*r0**2
c
c.... mass conversion (Mbox = rho0 * Lbox^3, Mbox_code = Ng^3
c     for non-cosmological run aM0 must be defined during initialization
c     [aM0] = [Msun]
      aM0 = rho0 * (boxh/hubble)**3 / ncell0
c
c.... constant L0 for the Sutherland & Dopita cooling function
c
c.... old
c      AL_SD = 2.0731d+22 * (4.d0 - 3.d0*Y_p) * (2.d0 - Y_p) * hubble /
c     &        sqrt(Om0) / r0**2  ! in ergs/s/cm^3
c.... new
c
       AL_SD = 1.6625d0 * (1.d0 - Y_p)**2 * hubble /
     &         sqrt(Om0) / r0**2  ! in ergs/s/cm^3
c
c.... constant for Compton cooling
c
      AL_Comp = 4.55d-4 * (4.d0 - 2.d0 * Y_p) / (8.d0 - 5.d0 * Y_p) * 
     &            T_CMB0**4 / hubble / sqrt(Om0)
c
#ifdef STARFORM 
      write (*,*)  'WARNING!'
      write (*,*)  'you may need some units for SF analysis'
      write (*,*)  ' E_0 , C_fb , C_SFR ,  rho_SF , in SF_Init () '
#endif
      write(*,*) 'set units to :'
      write (*,'(1x,''r0 [h^{-1} Mpc]    = '',1pe12.4)') r0
      write (*,'(1x,''t0 [yrs]           = '',1pe12.4)') t0
      write (*,'(1x,''rho0 [M_sun/Mpc^3] = '',1pe12.4)') rho0
      write (*,'(1x,''v0 [km/s]          = '',1pe12.4)') v0
      write (*,'(1x,''T0 [K]             = '',1pe12.4)') T_0
      write (*,'(1x,''P0 [g/cm/s^2]      = '',1pe12.4)') P0
      write (*,'(1x,''S0 [keV cm^2]      = '',1pe12.4)') S_0
      end
c-------------------------------------------------------------------- 
c     -------------------------
      subroutine Set_Cooling ()
c     -------------------------
c
c     tabulate cooling curve as a function of T in units of 10^4 K
c     table entry i corresponds to i = 100*T^1/4
c     AL_0 must be set in Set_Units prior to call to this routine
c      
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

#ifdef CLOUDY_COOLING
      include 'a_clcool.fh' 
      real*8 cdum, hdum, ct_crit 

#ifdef SD93_COOLING
c.... check that other cooling tables are not defined 
      call Open_ASCII_File ( iOErr , errorfile , sequent ) 
      write(iOErr,*) '1 : error : in Set_Cooling:'
      write(iOErr,*)
     &    ' both CLOUDY_COOLING and SD93_COOLING are defined in a_def.h'
      write(iOErr,*)
     &     ' while only one must be defined. stopping... '
      close (iOErr) 
      stop
#endif
c
c.... use net cooling/heating rate tabulated using Cloudy
c
      open ( 40 , file = 'clcool.dat' )
      read(40,*) tlmin, tlmax, dlt, nlt
      read(40,*) dlmin, dlmax, dld, nld
      read(40,*) Zlmin, Zlmax, dlZ, nlz
      read(40,*) rsmin, rsmax, drs, nrs

      write(*,*) tlmin, tlmax, dlt, nlt
      write(*,*) dlmin, dlmax, dld, nld
      write(*,*) Zlmin, Zlmax, dlZ, nlz
      write(*,*) rsmin, rsmax, drs, nrs

      dlti = 1.d0 / dlt
      dldi = 1.d0 / dld
      dlzi = 1.d0 / dlz
      drsi = 1.d0 / drs 

      do irs = 1, nrs
        do ilz = 1 , nlz
          do ild = 1 , nld 

            ct_crit = 0.0
            do ilt = 1 , nlt
              read(40,*) d1, d2, d3, d4, d5, d6, d7, d8, d9, cdum, hdum
              if ( d1 .ge. 3.2 .and. ct_crit .eq. 0.0 ) then
                ct_crit = cdum  
              endif
              cdum = max ( cdum , smallrate )
              cdum = max ( cdum , ct_crit   )  ! fix for a trough in equilibrium H2 cooling curve
              hdum = max ( hdum , smallrate ) 
              coolcl(ilt,ild,ilz,irs) = (cdum - hdum) * 1.d23 * AL_SD  ! COOLING-HEATING
c               coolcl(ilt,ild,ilz,irs) = (cdum ) * 1.d23 * AL_SD  ! COOLING ONLY
              f_ion(ilt,ild,ilz,irs) = d6 / 10.d0**(d2)
 100          format(6(e12.6,1x))
            enddo
          enddo
        enddo
      enddo
      close ( 40 ) 

      return
#endif

#ifdef SD93_COOLING
c
c.... use CIE cooling curves of Sutherland & Dopita 1993
c
      real csmall, Zsolar
      parameter ( ntab = 5000, nmetal = 8  ) 
      parameter ( csmall = 1.0e-10, Zsolar = 0.0199 )
      common / CLRATE / coolfun(0:ntab,nmetal), zval(nmetal)
      real xd, temp(200), coolr(200,nmetal), ap(nmetal)
      character*79 Line
      real CINTR
      real*8 T4, cr_H2, Z_met

      open (40 , file = 'coolfunc.dat' )
      read (40,40) Line
 40   format(A)
      write (*,*) 'setting cooling rates...'
      write (*,40) Line
      read (40,*) nmet, (zval(j),j=1,nmet)
      if (nmet .ne. nmetal) then 
        call Open_ASCII_File ( iOErr , errorfile , sequent ) 
        write(iOErr,*) '1 : *error in Set_Cooling: nmet .ne. nmetal'
        write(iOErr,*) 'nmet =',nmet,' nmetal =',nmetal
        write(iOErr,*) 'stopping...'
        close (iOErr)
        stop
      endif

      n = 0
 20   read(40,*,end=30,err=30) aT,(ap(j),j=1,nmetal)
      n = n + 1
      temp(n) =10.d0**(aT-4.d0)
      do j = 1 , nmetal
        coolr(n,j) = (10.d0**(ap(j)+23.0) * AL_SD)
      enddo
      go to 20

 30   write (*,*) ' number of lines read for cooling=', n
      do j = 1 , nmetal 
        Z_met = 10.d0**(zval(j))
        do i = 0 , ntab
          xd = (float(i)/100.0)**4
          coolfun(i,j)  = CINTR ( xd, temp, coolr(1,j), n )
c
c....     H2 cooling at low T
c          
          T4 = xd 
          call CoolingRate_H2 ( T4 , Z_met , cr_H2 )
          cr_H2 = cr_H2 * 1.d23 * AL_SD
          if ( cr_H2 .gt. coolfun(i,j) .and. T4 .lt. 4.25 ) then 
            coolfun(i,j) = cr_H2
          endif
        enddo 
      enddo
c
c.... small ad hoc cooling rate below 10^4K 
c
c      do j = 1 , nmetal
c        do i = 0 , 99
c          dee = min( ((float(i)-100.) * 0.2)**2, 30. )
c          coolfun(i,j)=coolfun(100,j)/exp(dee) +
c     &                 5.0e-10 * coolfun(100,j)
c        enddo
c      enddo
      
      close (40)
      return
#endif
      
      call Open_ASCII_File ( iOErr , errorfile , sequent ) 
      write(iOErr,*) '1 : error : in Set_Cooling:'
      write(iOErr,*)
     &     ' routine is called but cool. rates are not defined'
      write(iOErr,*)
     &     ' CLOUDY_COOLING or SD93_COOLING must be defined in a_def.h'
      close (iOErr) 
      stop
      
      end
c     ---------------------------------------------------
      subroutine CoolingRate_H2 ( T_g , Z_met, coolrate )
c     ---------------------------------------------------
c
c     approximation to cooling rate due to H2 
c     in the low-density regime 
c     [Galli & Palla 1998, A&A 335, 403 [eq. A7]]
c     in this regime (n_H <~ 10^2 cm^-3) cooling rate 
c     is independent of density. 
c        
c     T_g is gas temperature in units of 10^4 K      
c     Z_met - metallicity in solar units = [Fe/H] / [Fe/H]_solar
c     coolrate is in erg sm^3 s^-1
c
c     the rate assumes that for zero metallicity H2 reaches universal 
c     abundance of n_H2 / N_H = 1.e-3 to 1.e-1  [see Haiman and Oh 2002]
c     the assumed value is set in parameter H2_u below 
c
c     for non-zero metallicities the cooling rate is assumed to be boosted
c     with an hoc additive term which approximately mimicks 
c     results of Spaans & Norman 1997, ApJ 483, 87 
c     (see also Wada and Norman 2001, ApJ 547, 172)
c     
c     to get the cooling rate multiply coolrate by n_H^2
c
      real*8 T_g, Z_met, coolrate
      real*8 lt, lt2, tm
      real*8 H2_u  ! assumed universal fraction of H2
      parameter ( H2_u = 1.d-2 )
      real*8 sigt
      parameter ( sigt = 0.75 )

      tm = max(T_g, 1.3d-3) * 1.d4    ! no cooling below 13 Kelvin ...
      tm = min(tm, 1.d5)        ! fixes numerics
      lt = dlog10(tm)
      lt2 = lt * lt
c     low density limit from Galli and Palla (1998) 
c     multiplied by 1d23 
      coolrate = 
     &           H2_u * 
     &           10.d0**(-103.d0 + 97.59d0 * lt     - 
     &                           (48.05d0       -
     &                            10.8d0   * lt + 
     &                            0.9032d0 * lt2  ) * lt2 ) 
c
c.... add a boost due to metals to mimick Spaans cooling curve
c     [see Fig.1 in Wada & Norman 2001, ApJ 547, 172
c
      if ( T_g .lt. 5.d0 ) then 
        coolrate = coolrate + 
     &      2.5d-25* dexp(-5.d-1*(lt - 4.)**2/sigt**2) / sigt * 
     &      (1.d0 + 10.d0 * Z_met)
      endif

      return
      end      
c----------------------------------------------------------------
      subroutine Open_ASCII_File ( iunit , filename , accesstype )
c     ------------------------------------------------------------
c
c     filename is the string containing the name of the file to open
c              (should end with the space character)
c     accesstype is the string containing the type of file access 
c                (e.g., append)
c              (should end with the space character)
c
# include "a_def.h"
c
      integer iunit
      character*256 filename, accesstype
c
      nfn  = index(filename, ' ') - 1
      nat  = index(accesstype, ' ') - 1
c
c.... take into account different syntax on AIX 
c     (define or comment OS_AIX in a_def.h as needed)
c
#ifdef OS_AIX
      open ( iunit , file = filename(1:nfn), 
     &     position = accesstype(1:nat) )
#endif
#ifndef OS_AIX
      open ( iunit , file = filename(1:nfn), 
     &               access = accesstype(1:nat) )
#endif
c
      return
      end
c     ----------------------------------
      subroutine GetRow ( irow , ifile )
c     ----------------------------------
c
c     input:  irow - number of record
c             ifile - number of i/o unit (0-99)
c     nrecl - length of ROW block in words
c
      include 'a_tree.h'
      read  ( ifile , rec = irow ) recdat
      return
      end
C-------------------------------------------------------
C SUBROUTINE FOR SCANNING THE MESH
C-------------------------------------------------------
c     -----------------------------
      subroutine Get_MaxLevelNow ()
c     -----------------------------
c
      include 'a_tree.h'
      include 'a_control.h'

      MaxLevelNow = MinLevel
      if ( iNOLL(MinLevel+1) .eq. nil ) then 
        MaxLevelNow = MinLevel
        return
      else
        do Level = MaxLevel , MinLevel + 1 , -1
          if ( iNOLL(Level) .gt. nil ) then 
            MaxLevelNow = Level 
            return
          endif
        enddo
      endif
      MaxLevelNow = Level - 1
      return
      end
c     -----------------------------------------
      subroutine Select_Cells ( Level , nLevel ) 
c     -----------------------------------------
c
c     purpose: selects cells (every eights) of a specified Level
c
      include 'a_tree.h'
      include 'a_control.h'

      integer Level , nLevel 

      nLevel = iNOLL(Level)     ! get level boundaries for index array 

      icell  = iGet_Cell_Index ( iHOLL(Level) )
      do j = 1 , nLevel         
        iSelect(j) = icell
        iOct       = (icell - ncell0) / nchild + 1
        iCell      = (iOctLL1(iOct) - 1) * nchild + ncell0 + 1
      enddo
      return
      end
c     ---------------------------------
      function iGet_Cell_Index ( iOct ) 
c     ---------------------------------
c
c     returns index of the first cell in iOct
c
      include 'a_tree.h'
      integer iOct
      iGet_Cell_Index = (iOct - 1) * nchild + ncell0 + 1
      return
      end
c     -------------------------------------
      subroutine Ps ( iC , Posx,Posy,Posz )      ! finds position
c     -------------------------------------
c     purpose : returns coordinates of cell center
c     Input   : iC     - pointer to a cell 
c     Output  : Posxyz    - x,y,z positions
c
      include 'a_tree.h'
      integer iC
      integer idelta(8,3) 
      data      idelta / -1,  1, -1,  1, -1,  1, -1,  1,
     &                   -1, -1,  1,  1, -1, -1,  1,  1, 
     &                   -1, -1, -1, -1,  1,  1,  1,  1 /
      if ( iC .le. ncell0 ) then 
        iC_ = iC - 1 
        i = iC_ / ng2               ! ng2 = ng**2 
        j = ( iC_ - i*ng2 ) / ng 
        k =   iC_ - i*ng2 - j*ng
        Posx = float(i) + 1.5
        Posy = float(j) + 1.5
        Posz = float(k) + 1.5
      else
        iC_ = iC + nbshift
        iO = ishft ( iC_ , - ndim )
        id = ishft ( 1, MaxLevel - iOctLv(iO) )   
        j  = iC_ + 1 - ishft( iO , ndim )
        Posx   = d_x * (iOctPs(1,iO) + sign ( id , idelta(j,1) ))
        Posy   = d_x * (iOctPs(2,iO) + sign ( id , idelta(j,2) ))
        Posz   = d_x * (iOctPs(3,iO) + sign ( id , idelta(j,3) )) 
      endif
      return
      end
C---------------------
c	PROFILES
C---------------------
C
c     -------------------------------------------------------
      subroutine Find_Halo_Center ( xc, yc, zc, ifindcenter )
c     Finding a new center of the halo
c     -------------------------------------------------------
c
c     ifindcenter = 0 -- gas density
c                 = 1 -- DM density
c                 = 2 -- total density
c                 = 3 -- potential
c                >= 4 -- no new center
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc , yc , zc 
      integer ifindcenter
c

c      write(*,*) 'Find_Halo_Center...'
      call Get_MaxLevelNow ()
c
c.... find a new center 
c
      if ( ifindcenter .le. 3 ) then

        rhogmax = -1.e6  ! max. gas density 
        rhotmax = -1.e6  ! max. total density
        phimin  = 1.e6   ! minimum of potential
        rhodmmax = -1.e6 
        irhogmax = -1
        irhotmax = -1 
        iphimin  = -1 
        irhodmmax = -1
        
        rmax = 10.0/r0
  
        DO Level = MinLevel , MaxLevelNow  
          IF ( Level .eq. MinLevel ) THEN
          ELSE 
            call Select_Cells ( Level , nLevel ) 
            do ic1 = 1 , nLevel
              icell = iSelect(ic1)
              do ic2 = 0 , 7
                ic = icell + ic2
                if ( iOctCh(ic) .eq. nil ) then 
                  call Ps ( ic , Posx , Posy , Posz )
c                  call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz)
c                  rr = sqrt( (Posx-xc)**2 +
c       &                     (Posy-yc)**2 +
c       &                     (Posz-zc)**2 )
c                  if ( rr .lt. rmax ) then                   
                    dcell1 = var(1,ic)*2.0**(3.*Level) + 1.0 ! tot. density
                    dcell2 = var(2,ic)  ! potential 
                    dcell3 = hvar(1,ic) ! gas density
                    dcell4 = ref(ic) ! mass of DM + stars
c                    write(*,*) dcell1, dcell2, dcell3, dcell4
                    if ( dcell1 .gt. rhotmax ) then 
                      rhotmax  = dcell1 
                      irhotmax = ic 
                    endif
                    if ( dcell2 .lt. phimin ) then 
                      phimin  = dcell2
                      iphimin = ic 
                    endif
                    if ( dcell3 .gt. rhogmax ) then
                      rhogmax  = dcell3
                      irhogmax = ic
                    endif
                    if ( dcell4 .gt. rhodmmax ) then 
                      rhodmmax = dcell4
                      irhodmmax = ic
                    endif
c                  endif
                endif
              enddo
            enddo
          ENDIF
        ENDDO      
c       
        if ( ifindcenter .eq. 0 ) then 
          if ( irhogmax .eq. -1 ) then 
            write(*,*) 'error in Find_Halo_Center:'
            write(*,*) '   failed to find gas density maximum'
            write(*,*) '   will use input center '
          else
            call Ps ( irhogmax , xc1 , yc1 , zc1 )
            write(*,*) 'input center x, y, z:', xc, yc, zc 
            write(*,*) 'identified gas density maximum rhogmax=',rhogmax
            write(*,*) 'at x, y, z =', xc1, yc1, zc1        
            write(*,*) 'in the cell',  irhogmax    
            write(*,*) 'will use it to construct profiles'
            xc = xc1
            yc = yc1 
            zc = zc1
          endif
        endif
        if ( ifindcenter .eq. 1 ) then
          if ( irhodmmax .eq. -1 ) then
            write(*,*) 'error in Find_Halo_Center:'
            write(*,*) '   failed to find DM density maximum'
            write(*,*) '   will use input center '
          else
            call Ps ( irhodmmax , xc1 , yc1 , zc1 )
            write(*,*) 'input center x, y, z:', xc, yc, zc
            write(*,*)'identified DM density maximum rhodmmax=',rhodmmax
            write(*,*) 'at x, y, z =', xc1, yc1, zc1
            write(*,*) 'will use it to construct profiles'
            xc = xc1
            yc = yc1
            zc = zc1
          endif
        endif
        if ( ifindcenter .eq. 2 ) then  ! total density
          if ( irhotmax .eq. -1 ) then 
            write(*,*) 'error in Find_Halo_Center:'
            write(*,*) '   failed to find total density maximum'
            write(*,*) '   will use input center '
          else
            call Ps ( irhotmax , xc1 , yc1 , zc1 )
            write(*,*) 'input center x, y, z:', xc, yc, zc 
            write(*,*) 'identified tot density maximum rhotmax=',rhotmax
            write(*,*) 'at x, y, z =', xc1, yc1, zc1            
            write(*,*) 'will use it to construct profiles' 
            xc = xc1
            yc = yc1 
            zc = zc1
          endif
        endif
        if ( ifindcenter .eq. 3 ) then  ! potential 
          if ( iphimin .eq. -1 ) then 
            write(*,*) 'error in Find_Halo_Center:'
            write(*,*) '   failed to find potential minimum'
            write(*,*) '   will use input center '
          else
            call Ps ( iphimin , xc1 , yc1 , zc1 )
c            write(*,*) 'input center x, y, z:', xc, yc, zc 
c            write(*,*) 'identified minimum of the potential=',phimin
c            write(*,*) 'at x, y, z =', xc1, yc1, zc1            
c            write(*,*) 'will use it to construct profiles' 
            xc = xc1
            yc = yc1 
            zc = zc1
          endif
        endif
        
      endif

      return
      end
c   (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
c
c     -----------------
      FUNCTION RANDd(M)
c     -----------------
c     random number generator
c     initialize with arbitrary integer > 0
c     after that never change it
c     -----------------
      DATA LC,AM,KI,K1,K2,K3,K4,L1,L2,L3,L4
     +/453815927,2147483648.,2147483647,536870912,131072,256,16777216,
     +   4,16384,8388608,128/
      ML=M/K1*K1
      M1=(M-ML)*L1
      ML=M/K2*K2
      M2=(M-ML)*L2
      ML=M/K3*K3
      M3=(M-ML)*L3
      ML=M/K4*K4
      M4=(M-ML)*L4
      M5=KI-M
      IF(M1.GE.M5)M1=M1-KI-1
      ML=M+M1
      M5=KI-ML
      IF(M2.GE.M5)M2=M2-KI-1
      ML=ML+M2
      M5=KI-ML
      IF(M3.GE.M5)M3=M3-KI-1
      ML=ML+M3
      M5=KI-ML
      IF(M4.GE.M5)M4=M4-KI-1
      ML=ML+M4
      M5=KI-ML
      IF(LC.GE.M5)ML=ML-KI-1
      M=ML+LC
      RANDd=M/AM
      RETURN
      END
c     ---------------------------------------
      function iFindCell ( L , xp , yp , zp )
c     ---------------------------------------
c
c     purpose: finds a cell of Level or above enclosing 
c              point with coordinates xp, yp, zp
c     input  : Level, xp, yp, zp
c
      include 'a_tree.h'
      parameter ( ng1 = ng + 1 )
      integer iCell, L
      double precision xp, yp, zp, xx, yy, zz

      ip    = int(xp)
      jp    = int(yp)
      kp    = int(zp)
      if ( ip .eq. 0 ) ip = ng 
      if ( jp .eq. 0 ) jp = ng
      if ( kp .eq. 0 ) kp = ng
      if ( ip .eq. ng1 ) ip = 1
      if ( jp .eq. ng1 ) jp = 1
      if ( kp .eq. ng1 ) kp = 1

      iCell = iCellIndex ( ip , jp , kp )

      if ( L .eq. MinLevel ) then
        iFindCell = iCell
        return
      endif

      xx    = xp - ip
      yy    = yp - jp
      zz    = zp - kp

      iter  = 0
 1    iter  = iter + 1
      dxr   = ishft(1,iter)  ! = 2**iter
      iChild = 1 +     mod ( int ( xx * dxr ) , 2 )
     &           + 2 * mod ( int ( yy * dxr ) , 2 )
     &           + 4 * mod ( int ( zz * dxr ) , 2 )
      iCCh = iOctCh(iCell)
      if ( iCCh .gt. nil ) then 
        iCell = ncell0 + ishft ( iCCh - 1 , ndim ) + iChild 
      endif
      if ( (iter.lt.L) .and. (iOctCh(iCell).gt.nil) ) go to 1

      iFindCell = iCell
            
      return
      end
c     ------------------------------------
      function iCellIndex ( ic , jc , kc )  
c     ------------------------------------
c
c     returns index of parent cell from iPr 
c     which corresponds to given combination of i,j,k
c
c     (!) Note: algorithm of cell identification depends on the
c     order in which we assume cells are placed in iPr:
c     assumed order is  i,j,k. When k is changing most rapidly
c     (see make_root_mesh subroutine)
c
      include 'a_setup.h'

      iCellIndex = ((ic-1)*ng + (jc-1))*ng + kc

      return
      end
c     -----------------------------------------------------------
      subroutine Output_Average_Profiles ( xc , yc , zc , rvir, 
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, ifindcenter , 
     &                            lpointsout, l_force_center,
     &                            fprofname )
c     -----------------------------------------------------------
c
c     input:  
c
c      integer ibintype : 0 - logarithmic, 1 - even (now only log)
c                         
c      real    deltavir : virial overdensity (with respect to mean)
c      integer ifindcenter (see below) 
c      logical lpointsout  - if .true, output sample points in points.dat
c      logical l_force_center = T -> use input (xc,yc,zc) as center
c                               F -> determine the center using 
c                                 max. gas density (ifindcenter=0)   
c                                 max. DM density ( =1)
c                                 max. total density ( =2)
c                                 min. grav. potential (=3)
c     character*256 fprofname - full path and file name to output profile to 
c                               it must end with a space
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc , yc , zc , rmin , rmax , deltavir
      integer nrbin, ibintype, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 fprofname 
c
      real*8 a2b, b2a
      real a2t, anow, tnow
      integer m
      parameter ( npoints = 1000000 )
      real xpn(npoints), ypn(npoints), zpn(npoints)
      real pdm(npoints), pd(npoints), pp(npoints)
      real ptemp(npoints), pent(npoints)
      common / POINTS / xpn, ypn, zpn, pdm, pd, pp, ptemp, pent
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins),dgbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), dtibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins),vbin(0:nmaxbins)
      real ndmp(0:nmaxbins),entbin(0:nmaxbins)
      real tdum(8), entdum(8)
      integer nbin(nmaxbins)
      integer Level, icell, iPack(8)
      real vel(8), velx(8), vely(8), velz(8)
      double precision xx, yy, zz
      integer iIP(8,3)                     ! set of pointers to the leftmost
                                           ! cells - to be used as kernels
      data iIP   / 1 , 3 , 1 , 5 , 1 , 3 , 1 , 0 ,
     &             3 , 5 , 5 , 0 , 3 , 0 , 0 , 0 ,
     &             5 , 0 , 0 , 0 , 0 , 0 , 0 , 0  /
      integer nfpn

      nfpn = index ( fprofname , ' ' ) - 1

      rlmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlmin) / float(nrbin) 

      nbins = int((rlmax-rlmin)/drl) + 1

      a3 = aexpn**3
      a2 = aexpn**2
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 

      call Get_MaxLevelNow ()
c
c.... find a new center 
c
      if (  .not. l_force_center ) then 
         call Find_Halo_Center ( xc, yc, zc, ifindcenter )
      endif

      m = 28548535
      xmax = -1.e6
      ymax = xmax
      zmax = xmax
      xmin = 1.e6
      ymin = xmin
      zmin = xmin
c
c.... select random samples in the volume 
c
      do i = 1 , npoints
        rr    = rmax * RANDd ( m )
        theta = pi2 * RANDd ( m ) - pi
        phi   = pi2 * RANDd ( m )
        xpn(i) = rr * cos(theta) * cos(phi) + xc
        ypn(i) = rr * cos(theta) * sin(phi) + yc
        zpn(i) = rr * sin(theta) + zc
        if ( xpn(i) .lt. xmin ) xmin = xpn(i)
        if ( xpn(i) .gt. xmax ) xmax = xpn(i)
        if ( ypn(i) .lt. ymin ) ymin = ypn(i)
	if ( ypn(i) .gt. ymax ) ymax = ypn(i)
        if ( zpn(i) .lt. zmin ) zmin = zpn(i)
        if ( zpn(i) .gt. zmax ) zmax = zpn(i)
      enddo
c
      do i = 1 , npoints
        xx  = xpn(i)
        yy  = ypn(i)
        zz  = zpn(i)
        if ( xx .lt. 1.0 ) then 
          xx = xx + ng
        else 
          if ( xx .gt. narr ) then 
            xx = xx - ng
          endif
        endif
        if ( yy .lt. 1.0 ) then 
          yy = yy + ng
        else 
          if ( yy .gt. narr ) then 
            yy = yy - ng
          endif
        endif
        if ( zz .lt. 1.0 ) then 
          zz = zz + ng
        else 
          if ( zz .gt. narr ) then 
            zz = zz - ng
          endif
        endif
        
        idcell = iFindCell ( MaxLevel , xx , yy , zz )
c
c....   gas density
c
        pd(i) = hvar(1,idcell)
c
c....   pressure
c
        pp(i) = hvar(6,idcell)

        if ( hvar(1,idcell) .gt. 0. ) then
          T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
	  ptemp(i) = T_0 * T_code * a2i
          pent(i)  = S_0 * T_code / hvar(1,idcell)**gamma1
        else
          write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
        endif
      enddo
c
c.... construct profiles
c
      do i = 0 , nmaxbins
        rbin(i) = 0.
        dgbin(i) = 0.
        ddmbin(i) = 0.
        ddmibin(i) = 0.
        dtibin(i) = 0.
        pbin(i) = 0.
	tbin(i) = 0.
        entbin(i) = 0.
        nbin(i) = 0
        ndmp(i) = 0.
      enddo

      do i = 1 , np
        xx = x(i)
        yy = y(i)
        zz = z(i)
c....   correct coordinate to insure periodic bundaries        
        call PeriodicBoundaryD( xc, yc, zc, xx, yy, zz )
        rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          ddmbin(ibin) = ddmbin(ibin) + pw(i)
          ndmp(ibin) = ndmp(ibin) + 1
        endif
      enddo

      if ( lpointsout ) then 
        open ( 12 , file = 'points.dat' )
      endif
      do ii = 1 , npoints
        rp = sqrt((xpn(ii)-xc)**2+(ypn(ii)-yc)**2+(zpn(ii)-zc)**2)
        ibin = max ( int((log10(rp) - rlmin)/drl) + 1 , 0 )
        if ( ibin .le. nmaxbins ) then
          nbin(ibin)   = nbin(ibin) + 1
 	  rbin(ibin)   = rbin(ibin) + rp
          dgbin(ibin)  = dgbin(ibin) + pd(ii)
          pbin(ibin) = pbin(ibin) + pp(ii)
          tbin(ibin) = tbin(ibin) + ptemp(ii)
	  entbin(ibin) = entbin(ibin) + pent(ii)
          if ( lpointsout ) then 
            write(12,*) rp*r0, pd(ii), ptemp(ii),pent(ii)
          endif
        endif
      enddo
      if ( lpointsout ) then       
        close ( 12 )
      endif
c
c.... prepare profiles
c
      do i = 0 , nbins
	if ( nbin(i) .gt. 0 ) then
          rbin(i) = rbin(i) / nbin(i) * r0
          dgbin(i) = dgbin(i) / nbin(i)
c          ddmbin(i) = ddmbin(i) / nbin(i)
          pbin(i) = pbin(i) / nbin(i)
          tbin(i) = tbin(i)/nbin(i)
          entbin(i) = entbin(i)/nbin(i)
	endif
      enddo

      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.0
      dti = 0.0
      ivirflag = 0 

      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
        else
          rl = 10.**(rlmin + float(i-1)*drl)
        endif
        rr = 10.**(rlmin + float(i)*drl)
        volr = 4.18879 * (rr**3 - rl**3)        
	if ( volr .gt. 0. ) then
          dmi = dmi + ddmbin(i) 
          dti = dti + ddmbin(i) + dgbin(i)*volr 
          ddmbin(i) = ddmbin(i)/volr 
          ddmibin(i) = dmi / (4.18879 * rr**3)
          dtibin(i) = dti / (4.18879 * rr**3)
          if ( i .gt. 0 .and.
     &         dtibin(i) .lt. deltavir .and.
     &         dtibin(i-1) .ge. deltavir ) then 
            ivir = i
            rvir = (deltavir * (rr - rl) + rl*dtibin(i) - 
     &             rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
c            rvir = rr 
          endif
        else
	  write(*,*) 'bad volr =',volr,i, rl, rr
        endif
      enddo
c
      open ( 12 , file = fprofname(1:nfpn) )
c
      aMvir = deltavir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      write(12,10) xc , yc , zc  
      write(12,11) rmin , rmax , nrbin , ibintype 
      write(12,12) ivir, rvir*r0, aMvir, deltavir
      write(12,13) ifindcenter, l_force_center
      write(12,14)
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc) = ', 3(g13.7,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &        ' nrbin =',i5,' ibintype =',i1)
 12   format ('# ivir = ',i4,' Rvir =',f7.4,
     &        ' [h^-1 Mpc], Mvir =',g13.6,
     &        ' [h^-1 Msun] for Dvir(mean)=',f7.2)
 13   format ('# ifindcenter =',i2,' l_force_center=',l6)
 14   format ('# Columns:')
 15   format ('# <r> r_l r_m r_r delta_g delta_dm deltai_dm',
     &        ' deltai_tot  T_g  S_g     N_dm_p' )
 16   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
     &        '             K  keV cm^2 ')
      do i = 0 , nbins
        if ( i .eq. 0 ) then 
          rl = 0.0
          rmid = 0.5*10.**rlmin 
        else
          rl = 10.**(rlmin + float(i-1)*drl)
          rmid = 10.**(rlmin + float(i)*drl - 0.5*drl) * r0
        endif
        rr = 10.**(rlmin + float(i)*drl)
	if ( nbin(i) .gt. 0 ) then
          write(12,17)
     &         rbin(i), rl*r0, rmid, rr*r0, 
     &         dgbin(i),ddmbin(i),ddmibin(i), dtibin(i), 
     &         tbin(i),entbin(i),int(ndmp(i))
	endif
      enddo
 17   format(10(2x,g14.7),1x,i7)
      close ( 12 )
      write(*,*) 'done writing profiles. exiting...'
c
      return
      end
c     -----------------------------------------------------------
      subroutine PeriodicBoundary ( xc, yc, zc, Posx, Posy, Posz)
c     Take care of periodic boundary conditions for analysis
c     -----------------------------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc, yc, zc, Posx, Posy, Posz
c      double precision xc, yc, zc, Posx, Posy, Posz
c
      diff_x = xc - Posx
      diff_y = yc - Posy
      diff_z = zc - Posz
      corr_x = zero
      corr_y = zero
      corr_z = zero
      
      if ( abs(diff_x) .gt. nf67 ) then
        if ( diff_x .gt. 0. ) then
          corr_x = ng
        else
          corr_x = -ng
        endif
      endif
      if ( abs(diff_y) .gt. nf67 ) then
        if ( diff_y .gt. 0. ) then
          corr_y = ng
        else
          corr_y = -ng
        endif
      endif
      if ( abs(diff_z) .gt. nf67 ) then
        if ( diff_z .gt. 0. ) then
          corr_z = ng
        else
          corr_z = -ng
        endif
      endif
      
      Posx = Posx + corr_x
      Posy = Posy + corr_y
      Posz = Posz + corr_z

      return
      end
c
c     -----------------------------------------------------------
      subroutine PeriodicBoundaryD ( xc, yc, zc, xx, yy, zz)
c     Take care of periodic boundary conditions for analysis
c     -----------------------------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc, yc, zc
      double precision xx, yy, zz
c
      diff_x = xc - xx
      diff_y = yc - yy
      diff_z = zc - zz
      corr_x = zero
      corr_y = zero
      corr_z = zero
      
      if ( abs(diff_x) .gt. nf67 ) then
        if ( diff_x .gt. 0. ) then
          corr_x = ng
        else
          corr_x = -ng
        endif
      endif
      if ( abs(diff_y) .gt. nf67 ) then
        if ( diff_y .gt. 0. ) then
          corr_y = ng
        else
          corr_y = -ng
        endif
      endif
      if ( abs(diff_z) .gt. nf67 ) then
        if ( diff_z .gt. 0. ) then
          corr_z = ng
        else
          corr_z = -ng
        endif
      endif
      
      xx = xx + corr_x
      yy = yy + corr_y
      zz = zz + corr_z

      return
      end


c     ---------------------------------------------------------
      subroutine Cooling_Rate ( rhogl, T_g , Z_met , coolrate )
c     ---------------------------------------------------------
c     input:  
c         rhogl - log10(n_H), where n_H is hydrogen number density in cm^-3
c         T_g  - gas temperature in units of 10^4 K
c         Z_met - log10([rho_metals/rho_gas/solar)
c         if ENRICH is not defined Z_met will be ignored 
c         if SD93_COOLING is defined, rhogl will be ignored 
c
c     output: coolrate - net output cooling/heating rate in code units
c        de/dt = e - coolrate * rho_g**2 / aexpn
c        where e is internal energy per unit volume and rho_g is 
c        gas density (e and rho_g are in code units)
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      real*8 rhogl, rho_g, T_g, Z_met, coolrate

#ifdef CLOUDY_COOLING
      include 'a_clcool.fh' 
c
c.... use net cooling/heating rate tabulated using Cloudy
c
      real*8 Tlog
      real*8 ac, bc, ah, bh 

      Tlog = log10(T_g) + 4.d0 
      it1 = int((Tlog - tlmin)*dlti) + 1
      it2 = it1 + 1
      id1 = int((rhogl - dlmin)*dldi) + 1
      id2 = id1 + 1
#ifndef METALCOOLING 
      iz1 = 1
      iz2 = 1
#endif
#ifdef METALCOOLING
      iz1 = int((Z_met - zlmin)*dlzi) + 1
      iz2 = iz1 + 1
#endif
      it1 = max(it1,1) 
      it1 = min(it1,nlt)
      it2 = max(it2,1) 
      it2 = min(it2,nlt)
      id1 = max(id1,1) 
      id1 = min(id1,nld)
      id2 = max(id2,1) 
      id2 = min(id2,nld)
      iz1 = max(iz1,1) 
      iz1 = min(iz1,nlz)
      iz2 = max(iz2,1) 
      iz2 = min(iz2,nlz)

      td = tlmin + dlt * (it1 - 1) + dlt
      dd = dlmin + dld * (id1 - 1) + dld
      zd = zlmin + dlz * (iz1 - 1) + dlz
      t1 = (td - Tlog) * dlti
      d1 = 1.d0 - t1
      t2 = (dd - rhogl) * dldi
      d2 = 1.d0 - t2
      t3 = (zd - Z_met) * dlzi 
      d3 = 1.d0 - t3
c
c.... trilinear (CIC) interpolation 
c
      coolrate = t1*t2*t3 * ccl_rs(it1,id1,iz1) +
     &           d1*t2*t3 * ccl_rs(it2,id1,iz1) + 
     &           t1*d2*t3 * ccl_rs(it1,id2,iz1) + 
     &           d1*d2*t3 * ccl_rs(it2,id2,iz1) +
     &           t1*t2*d3 * ccl_rs(it1,id1,iz2) + 
     &           d1*t2*d3 * ccl_rs(it2,id1,iz2) + 
     &           t1*d2*d3 * ccl_rs(it1,id2,iz2) + 
     &           d1*d2*d3 * ccl_rs(it2,id2,iz2)

      return
#endif

#ifdef SD93_COOLING
c
c.... use CIE cooling curves of Sutherland & Dopita 1993
c
      real*8 dummy, Z1
      real*8 csmall, Zsolar, Zsolari, Zpri, slope, zerop
      parameter ( ntab = 5000 , nmetal = 8 , csmall=1.d-10 )
      parameter ( Zsolar = 0.0199 , Zsolari = 1.d0/Zsolar )
      parameter ( Zpri = 1.d-5 * Zsolar )
      common / CLRATE / coolfun(0:ntab,nmetal), zval(nmetal)

      itab = int(100.0 * sqrt(sqrt(max(T_g,1.d-4))))
      itab = max( itab , 0 )
      itab = min( itab , ntab )
#ifndef METALCOOLING   
      coolrate = coolfun(itab,1) ! primordial composition 
#endif
#ifdef METALCOOLING 
c
c.... in the following, the small difference between primordial and Z/Zsun=-3
c     will be ignored 
c
      if( Z_met .gt. Zval(2) ) then
         Z1 = Z_met
      else
        coolrate = coolfun(itab,1) ! primordial composition
        return
      endif
c
c...  if Z > Fe/H=0.5, extrapolate linearly in log10(Z/Zsun)
c
      if ( Z1 .gt. zval(nmetal) ) then
        slope =  (coolfun(itab,nmetal) - coolfun(itab,nmetal-1)) /
     &           (zval(nmetal) - zval(nmetal-1))
        zerop = coolfun(itab,nmetal) - slope * zval(nmetal)
        coolrate = slope * Z1 + zerop
        return
      endif
c
c.... interpolate linearly in log(Z/Zsun)
c
      do j = 3 , nmetal
        if ( Z1 .le. zval(j) ) then
          slope = (coolfun(itab,j)-coolfun(itab,j-1)) / 
     &            (zval(j)-zval(j-1))
          zerop = coolfun(itab,j) - slope*zval(j)
          coolrate = slope * Z1 + zerop
          return
        endif
      enddo

      call Open_ASCII_File ( iOErr , errorfile , sequent ) 
      write(iOErr,*) '1 : error : in Cooling_Rate: something is wrong:'
      write(iOErr,*) 'T_g, Z_met =', T_g, Z_met
      close (iOErr) 
      stop
#endif
#endif

      return
      end
c     -------------------------------------
      subroutine Set_Cooling_Rate_rs ( rs ) 
c     -------------------------------------
c
c     prepare cooling rate table for a given redshift rs
c     from the Cloudy cooling rate table 
c
c     real*8 rs - redshift 
c     
c     this routine is only used if CLOUDY_COOLING is defined in a_def.h
c     it should be called only in the beginning of every step 
c
c
      include 'a_clcool.fh' 
      real*8 rs, rsd
      real*8 ac, bc, ah, bh

      irs = int((rs - rsmin)*drsi) + 1
      irs1 = max(irs,1)
      irs1 = min(irs1,nrs)
      irs2 = min(irs+1,nrs)
      irs2 = max(irs+1,1)

      if ( irs1 .eq. irs2 ) then 
        do ilt = 1 , nlt
        do ild = 1 , nld
        do ilz = 1 , nlz 
          ccl_rs(ilt,ild,ilz) = coolcl(ilt,ild,ilz,irs1)
        enddo
        enddo
        enddo
      else
        rs1 = rsmin + drs*(irs1-1)
        rs2 = rsmin + drs*(irs2-1)

        do ilt = 1 , nlt
        do ild = 1 , nld
        do ilz = 1 , nlz 
          ac = (coolcl(ilt,ild,ilz,irs2) - coolcl(ilt,ild,ilz,irs1)) /
     &         (rs2 - rs1)
          bc = coolcl(ilt,ild,ilz,irs1) - ac * rs1 
          ccl_rs(ilt,ild,ilz) = ac * rs + bc 
        enddo
        enddo
        enddo
      endif

      return
      end
c-------------------------------------------------------------------- 
c
c     -------------------------------------
      FUNCTION zbrent(func,fp,np,x1,x2,tol)
c     -------------------------------------
c
      INTEGER ITMAX
      REAL zbrent
      integer np 
      real*8 fp(np)
      real*8 tol,x1,x2,EPS
      real*8 func
      EXTERNAL func
      PARAMETER (ITMAX=100,EPS=3.e-8)
      INTEGER iter
      REAL*8 a,b
      REAL c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      a=x1
      b=x2
      fa=func(a,fp,np)
      fb=func(b,fp,np)
      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then 
        write(*,*) 'root must be bracketed for zbrent'
        write(*,*) 'fa =',fa,'  fb =', fb
        pause
      endif
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc).lt.abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1 .or. fb.eq.0.)then
          zbrent=b
          return
        endif
        if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
          s=fb/fa
          if(a.eq.c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p.gt.0.) q=-q
          p=abs(p)
          if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d) .gt. tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=func(b,fp,np)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END
c$$$C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
c$$$c     --------------------------------------------------------------
c$$$      double precision function INTEGRATE(FUNC,fp,np,a,b,dxinit,eps)
c$$$c     --------------------------------------------------------------
c$$$c
c$$$c     Quadrature using fifth order Runge-Kutta with adaptive step size.
c$$$c     Based on Press et al, Numerical Recipes in C, 2nd ed, pp 719-722.
c$$$c
c$$$c     Runge-Kutta driver with adaptive stepsize control.  Integrate starting
c$$$c     value y from a to b with accuracy eps, storing intermediate results in
c$$$c     global variables.  dxinit should be set as a guessed first stepsize.
c$$$c
c$$$c     Pass a vector of parameters of length np to FUNC in fp(np).
c$$$c
c$$$c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c$$$c     personal, educational or research purposes.  Commercial purposes require
c$$$c     special arrangements. If you publish a paper that depends upon this code,
c$$$c     please cite it appropriately.
c$$$c
c$$$c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c$$$c
c$$$c     slight modification by A.Kravtsov to allow to pass a vector 
c$$$c     of parameters to FUNC
c$$$
c$$$      implicit none
c$$$      integer np 
c$$$      double precision a, b, eps, dxinit, FUNC, fp(np)
c$$$      external FUNC
c$$$
c$$$      integer maxsteps
c$$$      parameter(maxsteps=100000000)
c$$$
c$$$      double precision x, dx, dxnext, y, dydx, yscale
c$$$      integer  Nstep
c$$$
c$$$      x     = a
c$$$      dx    = dxinit
c$$$      y     = 0.d0
c$$$      Nstep = 0
c$$$
c$$$      do while ((x-b)*(b-a).lt.0.d0.and.Nstep.lt.maxsteps)
c$$$        Nstep = Nstep + 1
c$$$        dydx = FUNC(x,fp,np)
c$$$c
c$$$c       yscale is the scaling used to monitor accuracy.  This general-purpose
c$$$c       choice can be modified if need be.
c$$$c
c$$$        yscale = max(abs(y) + abs(dx*dydx), 1.d-20)
c$$$        if ((x+dx-b)*(x+dx-a).gt.0.d0)  ! If stepsize overshoots, decrease it.
c$$$     1    dx = b - x
c$$$
c$$$        call RUNGE5VAR(y,dydx,x,dx,eps,yscale,dxnext,FUNC,fp,np)
c$$$
c$$$        dx = dxnext
c$$$      end do
c$$$
c$$$      if (Nstep.ge.maxsteps)
c$$$     1  write (*,*) 'WARNING: failed to converge in INTEGRATE.'
c$$$
c$$$      INTEGRATE = y
c$$$
c$$$      return
c$$$      end
c$$$c
c$$$c     -------------------------------------------------------------
c$$$      SUBROUTINE RUNGE5VAR(y,dydx,x,htry,eps,yscale,hnext,DERIVS,
c$$$     1                     fp,np)
c$$$c     -------------------------------------------------------------
c$$$c
c$$$c     Fifth-order Runge-Kutta step with monitoring of local truncation error
c$$$c     to ensure accuracy and adjust stepsize.  Input are the dependent
c$$$c     variable y and its derivative dydx at the starting value of the
c$$$c     independent variable x.  Also input are the stepsize to be attempted
c$$$c     htry, the required accuracy eps, and the value yscale, against which the
c$$$c     error is scaled.  On output, y and x are replaced by their new values.
c$$$c     hdid is the stepsize that was actually accomplished, and hnext is the
c$$$c     estimated next stepsize.  DERIVS is the user-supplied routine that
c$$$c     computes right-hand-side derivatives.  The argument fp is a vector 
c$$$c     of parameters (np parameters) to be passed to DERIVS 
c$$$c     (NOT integrated over).
c$$$c
c$$$c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c$$$c     personal, educational or research purposes.  Commercial purposes require
c$$$c     special arrangements. If you publish a paper that depends upon this code,
c$$$c     please cite it appropriately.
c$$$c
c$$$c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c$$$c
c$$$c     slight modification by A.Kravtsov to allow to pass a vector 
c$$$c     of parameters to FUNC
c$$$c
c$$$      implicit none
c$$$      integer np 
c$$$      double precision eps,hnext,htry,x,dydx,y,yscale,DERIVS,fp(np)
c$$$      external DERIVS
c$$$
c$$$      double precision errmax,h,hold,htemp,xnew,yerr,ytemp
c$$$
c$$$      double precision safety,pgrow,pshrink,errcon
c$$$      parameter (safety  =  0.9d0)
c$$$      parameter (pgrow   = -0.2d0)
c$$$      parameter (pshrink = -0.25d0)
c$$$      parameter (errcon  =  1.89d-4)
c$$$
c$$$      h = htry                         ! Set stepsize to initial accuracy.
c$$$      errmax = 10.d0
c$$$      do while (errmax.gt.1.d0)
c$$$        call RUNGE(y,dydx,x,h,ytemp,yerr,DERIVS,fp,np)
c$$$
c$$$        errmax = abs(yerr/yscale)/eps   ! Scale relative to required accuracy.
c$$$        if (errmax.gt.1.d0) then        ! Truncation error too large; reduce h
c$$$          htemp = safety*h*(errmax**pshrink)
c$$$          hold = h
c$$$          h = sign(max(abs(htemp),0.1d0*abs(h)),h)  ! No more than factor of 10
c$$$          xnew = x + h
c$$$          if (xnew.eq.x) then
c$$$            write (*,*) 'WARNING: ',
c$$$     1                  'Stepsize underflow in RUNGE5VAR().'
c$$$            h = hold
c$$$            errmax = 0.d0
c$$$          end if
c$$$        end if
c$$$      end do
c$$$c
c$$$c     Step succeeded.  Compute estimated size of next step.
c$$$c
c$$$      if (errmax.gt.errcon) then
c$$$        hnext = safety*h*(errmax**pgrow)
c$$$      else
c$$$        hnext = 5.d0 * h                ! No more than factor of 5 increase.
c$$$      end if
c$$$      x = x + h
c$$$
c$$$      y = ytemp
c$$$
c$$$      return
c$$$      end
c$$$c
c$$$c     ---------------------------------------------------
c$$$      SUBROUTINE RUNGE(y,dydx,x,h,yout,yerr,DERIVS,fp,np)
c$$$c     ---------------------------------------------------
c$$$c
c$$$c     Given values for a variable y and its derivative dydx known at x, use
c$$$c     the fifth-order Cash-Karp Runge-Kutta method to advance the solution
c$$$c     over an interval h and return the incremented variables as yout.  Also
c$$$c     return an estimate of the local truncation error in yout using the
c$$$c     embedded fourth order method.  The user supplies the routine
c$$$c     DERIVS(x,y,dydx), which returns derivatives dydx at x.
c$$$c
c$$$c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c$$$c     personal, educational or research purposes.  Commercial purposes require
c$$$c     special arrangements. If you publish a paper that depends upon this code,
c$$$c     please cite it appropriately.
c$$$c
c$$$c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c$$$c
c$$$c     slight modification by A.Kravtsov to allow to pass a vector 
c$$$c     of parameters to FUNC
c$$$c
c$$$      implicit none
c$$$
c$$$      integer np
c$$$      double precision h,x,dydx,y,yerr,yout,DERIVS,fp(np)
c$$$
c$$$      external DERIVS
c$$$
c$$$      double precision ak3, ak4, ak5 ,ak6
c$$$
c$$$      double precision a2,a3,a4,a5,a6
c$$$      double precision c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
c$$$      parameter(a2  =    0.2d0)
c$$$      parameter(a3  =    0.3d0)
c$$$      parameter(a4  =    0.6d0)
c$$$      parameter(a5  =    1.d0)
c$$$      parameter(a6  =    0.875d0)
c$$$      parameter(c1  =   37.d0/378.d0)
c$$$      parameter(c3  =  250.d0/621.d0)
c$$$      parameter(c4  =  125.d0/594.d0)
c$$$      parameter(c6  =  512.d0/1771.d0)
c$$$      parameter(dc1 = c1 -  2825.d0/27648.d0)
c$$$      parameter(dc3 = c3 - 18575.d0/48384.d0)
c$$$      parameter(dc4 = c4 - 13525.d0/55296.d0)
c$$$      parameter(dc5 = -277.d0/14336.d0)
c$$$      parameter(dc6 = c6 -     0.25d0)
c$$$
c$$$      ak3 = DERIVS(x+a3*h,fp,np)
c$$$      ak4 = DERIVS(x+a4*h,fp,np)
c$$$      ak5 = DERIVS(x+a5*h,fp,np)
c$$$      ak6 = DERIVS(x+a6*h,fp,np)
c$$$c
c$$$c     Estimate the fifth order value.
c$$$c
c$$$      yout = y + h*(c1*dydx + c3*ak3 + c4*ak4  + c6*ak6)
c$$$c
c$$$c     Estimate error as difference between fourth and fifth order
c$$$c
c$$$      yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)
c$$$
c$$$      return
c$$$      end
c     --------------------------------------------------
      subroutine X_ion ( rhogl , T_g , Z_met , rs , xion )
c     --------------------------------------------------
c
c     equilibrium ionized fraction for a given condition
c     calculated using Cloudy tables
c
c     input:  
c         rhogl - log10(n_H), where n_H is hydrogen number density in cm^-3
c         T_g  - gas temperature in K
c         Z_met - log10([rho_metals/rho_gas/solar)
c         if ENRICH is not defined Z_met will be ignored 
c
c     output: xion - ionized fraction = n_HII/n_H
c
#     include "a_def.h"
c
      include 'a_clcool.fh' 
      include 'a_tree.h'
      include 'a_control.h'

      real*8 rs, rhogl, T_g, Z_met, xion 
      real*8 ac, bc, xion1, xion2

      Tlog = log10(T_g)  
      it1 = int((Tlog - tlmin)*dlti) + 1
      it2 = it1 + 1
      id1 = int((rhogl - dlmin)*dldi) + 1
      id2 = id1 + 1
#ifndef METALCOOLING 
      iz1 = 1
      iz2 = 1
#endif
#ifdef METALCOOLING
      iz1 = int((Z_met - zlmin)*dlzi) + 1
      iz2 = iz1 + 1
#endif
      it1 = max(it1,1) 
      it1 = min(it1,nlt)
      it2 = max(it2,1) 
      it2 = min(it2,nlt)
      id1 = max(id1,1) 
      id1 = min(id1,nld)
      id2 = max(id2,1) 
      id2 = min(id2,nld)
      iz1 = max(iz1,1) 
      iz1 = min(iz1,nlz)
      iz2 = max(iz2,1) 
      iz2 = min(iz2,nlz)

      td = tlmin + dlt * (it1 - 1) + dlt
      dd = dlmin + dld * (id1 - 1) + dld
      zd = zlmin + dlz * (iz1 - 1) + dlz
      t1 = (td - Tlog) * dlti
      d1 = 1.d0 - t1
      t2 = (dd - rhogl) * dldi
      d2 = 1.d0 - t2
      t3 = (zd - Z_met) * dlzi 
      d3 = 1.d0 - t3

      irs = int((rs - rsmin)*drsi) + 1
      irs1 = max(irs,1)
      irs1 = min(irs1,nrs)
      irs2 = min(irs+1,nrs)
      irs2 = max(irs+1,1)

      xion1 = t1*t2*t3 * f_ion(it1,id1,iz1,irs1) +
     &        d1*t2*t3 * f_ion(it2,id1,iz1,irs1) + 
     &        t1*d2*t3 * f_ion(it1,id2,iz1,irs1) + 
     &        d1*d2*t3 * f_ion(it2,id2,iz1,irs1) +
     &        t1*t2*d3 * f_ion(it1,id1,iz2,irs1) + 
     &        d1*t2*d3 * f_ion(it2,id1,iz2,irs1) + 
     &        t1*d2*d3 * f_ion(it1,id2,iz2,irs1) + 
     &        d1*d2*d3 * f_ion(it2,id2,iz2,irs1)

      xion2 = t1*t2*t3 * f_ion(it1,id1,iz1,irs2) +
     &        d1*t2*t3 * f_ion(it2,id1,iz1,irs2) + 
     &        t1*d2*t3 * f_ion(it1,id2,iz1,irs2) + 
     &        d1*d2*t3 * f_ion(it2,id2,iz1,irs2) +
     &        t1*t2*d3 * f_ion(it1,id1,iz2,irs2) + 
     &        d1*t2*d3 * f_ion(it2,id1,iz2,irs2) + 
     &        t1*d2*d3 * f_ion(it1,id2,iz2,irs2) + 
     &        d1*d2*d3 * f_ion(it2,id2,iz2,irs2)


      if ( irs1 .eq. irs2 ) then 
         xion = xion1
      else
        rs1 = rsmin + drs*(irs1-1)
        rs2 = rsmin + drs*(irs2-1)
        ac = (xion2 - xion1) /
     &         (rs2 - rs1)
        bc = xion1 - ac * rs1 
        xion = ac * rs + bc 
      endif
      return
      end

c     --------------------------------------------
      subroutine Read_Gas_Binary_HEADER ( dt , FileName )
c     --------------------------------------------
c
c     purpose: reads files with control information (no data)

      include 'a_tree.h'
      include 'a_control.h'
      real*8 dt
      character*256 FileName
      character*256 jname
      real*8 b2a

      nfname = index ( FileName , ' ' ) - 1
      open  ( 19 , file = FileName(1:nfname) , form = 'unformatted' )
      write(*,*) 'reading tree data from file ',FileName(1:nfname)
      read ( 19 ) jname 
      if ( jname .ne. jobname1 ) then 
        write(*,*) '* warning : Read_Gas_Binary : jname .ne. jobname1'
        write(*,'(2(A256))') jname, jobname1
c        stop
      endif
      read ( 19 ) istep , t , dt, aexpn, ainit  
      write(*,20) 
     &      istep, t, dt, aexpn, ainit
 20   format ('istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)
      do Level = MinLevel , MaxLevel 
        aexp(Level) = aexpn
      enddo
c     if having trouble on linux, try inserting a dummy first in the read
c     list below
      read ( 19 )  boxh, Om0, Oml0, Omb0, hubble
      write(*,21) 
     &      boxh, Om0, Oml0, Omb0, hubble      
 21   format('Lbox[/h Mpc] =',f9.4,1x,' Om0 = ',f5.3,1x,
     &        ' Oml0= ',f5.3,1x,' Omb0= ',f7.5,1x,' h = ',f5.3)
      read ( 19 ) nextras
      if ( nextras .gt. nextra ) then 
        write(*,*) '* error: Read_Gas_Binary : nextras > nextra :',
     &             nextras, nextra
        write(*,*) '* (check nextra parameter in a_control.h) '
        stop
      endif
      read ( 19 ) (extra(i), i=1,nextras)
      read ( 19 ) (lextra(i), i=1,nextras)

      read ( 19 ) MinLev , MaxLevelNow
      if ( MaxLevelNow .gt. MaxLevel ) then 
         write(*,*)
     &        '* error: Read_Gas_Binary : MaxLevelNow > MaxLevel :',
     &        MaxLevelNow, MaxLevel 
         stop
      endif
      close (19)
      write(*,*) 'done reading header of hydro file'
      end
c     -------------------
      function iLv ( iC )      ! 
c     -------------------
c
c     purpose : finds level of cell iC
c     Input   : iC     - pointers to cells (iC<=ncell0 --> iLv=0)

      include 'a_tree.h'
      integer iC
      if ( iC .le. ncell0 ) then 
        iLv = MinLevel
      else
        iLv = iOctLv(ishft ( iC + nbshift , - ndim ))
      endif
      return
      end
c     ---------------------------
      real*8 function age ( td )
c     ---------------------------
c
c     returns age of the Universe at td (time in code units)
c
c     uses Oleg's formula for flat LCDM 
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 td, b2a

      as = b2a ( td ) ! convert code time to expansion factor
      ff = Om0/(1.-Om0)/as**3 
c
c  .... calculate age of the universe in Gyrs 
c
      age = 9.779/hubble * 2./3./
     &      sqrt(1.-Om0)*log((1.+sqrt(1.+ff))/sqrt(ff))

      return
      end
