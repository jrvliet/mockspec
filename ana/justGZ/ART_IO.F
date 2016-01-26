c     =============================================================
c
c      ART Version 4: ART_IO.f - routines dealing with Input/Output
c
c     =============================================================
c     ------------------------------------------------------------
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
#ifdef OS_AIX
      open ( iunit , file = TRIM(filename), 
     &     position = TRIM(accesstype) )
#endif
#ifndef OS_AIX
      open ( iunit , file = TRIM(filename), 
     &               access = TRIM(accesstype) )
#endif
c
      return
      end
c
c     --------------------------------------------
      subroutine Read_Gas_Binary ( dt , FileName )
c     --------------------------------------------
c
c     purpose: reads files with control information and gas data
c              and hydrodynamic variables for cells
c
      include 'a_tree.h'
      include 'a_control.h'
      real*8 dt
      character*256 FileName
      character*256 jname
      real*8 b2a
c
      open (iOIO, file = TRIM(iofile), position='append')
c      
      open  ( 19 , file = TRIM(FileName) , form = 'unformatted' )
      write(iOIO,*) 'reading tree data from file ',TRIM(FileName)
      read ( 19 ) jname 
      if ( jname .ne. jobname1 ) then 
        write(iOIO,*) '* warning: Read_Gas_Binary : jname .ne. jobname1'
        write(iOIO,'(2(A256))') jname, jobname1
c        stop
      endif
      read ( 19 ) istep , t , dt, aexpn, ainit  
      write(iOIO,20) 
     &      istep, t, dt, aexpn, ainit
 20   format ('istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)
      do Level = MinLevel , MaxLevel 
        aexp(Level) = aexpn
      enddo
c     if having trouble on linux, try inserting a dummy first in the read
c     list below
      read ( 19 ) boxh, Om0, Oml0, Omb0, hubble
      write(iOIO,21) 
     &      boxh, Om0, Oml0, Omb0, hubble      
 21   format('Lbox[/h Mpc] =',f9.4,1x,' Om0 = ',f6.3,1x,
     &        ' Oml0= ',f6.3,1x,' Omb0= ',f8.5,1x,' h = ',f6.3)
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
        write(iOIO,*) 'iOctFree =',iOctFree,' nOct =', nOct
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
            write(*,*) 'error : L =',Lev,' nLevCells =',nLevCells,
     &        ' > nclmax =',nclmax,' set in a_setup.h'
            write(*,*) '=> increase nclmax and rerun.'
            close ( 19 )
            stop
          endif
          ntot = ntot + nLevel 
          !write(iOIO,*) 'reading tree (oct) data for level ', Lev
          do ic1 = 1 , nLevel
            read(19) (iOctPs(i,iOct),i=1,3),(iOctNb(i,iOct),i=1,6),
     &                iOctPr(iOct), iOctLv(iOct), iOctLL1(iOct), 
     &                iOctLL2(iOct)
            iOct = iOctLL1(iOct)
          enddo
          write(iOIO,*) 'reading cell data: ncells =', nLevel*nchild
          do ic1 = 1 , nLevel*nchild
            read ( 19 ) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar), 
     &                 (var(i,idc), i=2,3)
          enddo
        enddo
      endif
      close (19)
      write(iOIO,*) 'done reading',ncell0+nchild*ntot,' cells...'

c     set previous time moment for the energy cons. control
      ap0    = b2a ( (1.d0*(t-5.d-1*dt)) )      
c
      close (iOIO)
      return
      end
c
c     ------------------------------------------------------
      subroutine Read_Gas_Binary_Adiabatic ( dt , FileName )
c     ------------------------------------------------------
c
c     purpose: reads files with control information and gas data
c              and hydrodynamic variables for cells
c
c     Only difference from Read_Gas_Binary () is 
c         nhvar - > nhvar-nchem
c
      include 'a_tree.h'
      include 'a_control.h'
      real*8 dt
c     Don't forget to change this to real*8
      character*256 FileName
      character*256 jname
      real*8 b2a
c
      open  ( 19 , file = TRIM(FileName) , form = 'unformatted' )
      write(*,*) 'reading tree data from file ',FileName(1:nfname)
      read ( 19 ) jname 
      if ( jname .ne. jobname1 ) then 
        write(*,*) '* warning : Read_Gas_Binary : jname .ne. jobname1'
        write(*,'(2(A256))') jname, jobname1
c        stop
      endif
      read ( 19 ) istep , t , dt, aexpn, ainit  
c     quick fixes to match types 
c      t = t1
c      dt = dt1
      write(*,20) 
     &      istep, t, dt, aexpn, ainit
 20   format ('istep =',i4,1x,' t =',g15.8,1x,
     &        'dt =',g15.8,1x,' a =',g15.8,' ainit =',g15.8)
      do Level = MinLevel , MaxLevel 
        aexp(Level) = aexpn
      enddo
c     if having trouble on linux, try inserting a dummy first in the read
c     list below
      read ( 19 ) boxh, Om0, Oml0, Omb0, hubble
      write(*,21) 
     &      boxh, Om0, Oml0, Omb0, hubble      
 21   format('Lbox[/h Mpc] =',f9.4,1x,' Om0 = ',f6.3,1x,
     &        ' Oml0= ',f6.3,1x,' Omb0= ',f8.5,1x,' h = ',f6.3)
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
      read ( 19 ) ((hvar(ivar,i),ivar=1,nhvar-nchem), i=1,ncell)
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
            write(*,*) 'error : L =',Lev,' nLevCells =',nLevCells,
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
            read ( 19 ) idc, iOctCh(idc), (hvar(i,idc),i=1,nhvar-nchem), 
     &                 (var(i,idc), i=2,3)
          enddo
        enddo
      endif
      close (19)
      write(*,*) 'done reading',ncell0+nchild*ntot,' cells...'

c     set previous time moment for the energy cons. control
      ap0    = b2a ( (1.d0*(t-5.d-1*dt)) )      
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
c
      open (iOIO, file = TRIM(iofile), position='append')
c
      open ( 3 , file = TRIM(FileName1) , form = 'unformatted' )

c.... read control information and check whether it has proper structure

      read      (3) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,hubble,Wp5,
     &                  Ocurv,Omb0,extras
      write (iOIO,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRIDC,NRECL,Om0,Oml0,Omb0,hubble
 100  format (1X,'Header=>',A45,/
     &           1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &           1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &           1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I6,/
     &           1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,
     &               ' Omegab_0=',F7.3,' Hubble=',f7.3)
      
c...  Working with new ICs (Daisuke added a line below)
      write(iOIO,*) ' nspecies =',nspecies
      write(iOIO,*) ' nspec    =',nspec
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

      open ( 1 , file = TRIM(FileName2), access = 'direct',
     &	         status = 'unknown', recl = nacces      )
 
      rewind 3


      N_particles = lspecies(nspecies)
      np = N_particles 

      write(iOIO,*) ' N Particles :',N_particles

#ifdef WORKOUT
      call Open_ASCII_File ( iOWork , workfile , append ) 
      do i = 1 , nspecies
        write(iOWork,*) 'ispecie =',i,lspecies(i)
      enddo
      write(iOWork,*) 'Read_Particles_Binary: reading np=',np,
     &                ' particles'
      close(iOWork)
#endif

      if ( np .eq. 0 .or. np .gt. npmax ) then 
        call Open_ASCII_File ( iOErr , errorfile , append ) 
        write(iOErr,*) '1 : error in Read_Particles_Binary: bad np =',
     &             np
        write(iOErr,*) 'stopping...'
        close ( iOErr )
        stop
      endif

      Npages      = (N_particles -1)/npage + 1
      N_in_last   = N_particles - npage*(Npages-1)
      write (iOIO,*) ' Pages=',Npages,' Species=',nspecies
      write (iOIO,*) ' NrowR=',nrowreal,' Nrow=',nrow
      write (iOIO,*) ' N_in_last=',N_in_last

      do irow = 1 , Npages         ! loop over particle pages
        In_page = npage
        if ( irow .eq. Npages ) In_page = N_in_last
         write (iOIO,*)' Read page=',IROW,' file=',ifile,' N=',In_page
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

C$OMP PARALLEL DO DEFAULT(SHARED)
      do ic1 = 1 , np 
        pt(ic1) = t   ! set t to the current time moment
      enddo
cc
cc.... compute each particle's weight
cc
c      n_spec = nspecies
c #ifdef STARFORM
c      n_spec = nspecies-1
c #endif
c      do ic0 = 1 , n_spec 
c        do ic1 = nsp(ic0,1) , nsp(ic0,2)
c           pw(ic1) = wpar(ic0)
c        enddo
c      enddo 

#ifdef STARFORM
      if ( start ) then         
        if ( nspecies+1 .gt. nspec ) then
          write (*,*)
     &        ' nspec in PARAMETER < nspecies+1 for starformation:'
          write (*,*) ' nspecies+1=',nspecies+1,' > nspec =',nspec
          STOP
        endif
        nspecies = nspecies + 1
        nsp(nspecies,1) = lsp(nspecies-1) + 1
        nsp(nspecies,2) = lsp(nspecies-1)
        lsp(nspecies) = nsp(nspecies,2)
        lspecies(nspecies) = lsp(nspecies)


        write(iOIO,*) 'in reading particles:'
        write(iOIO,*) '    lsp      wpar         nsp1      nsp2 '
        do i = 1 , nspecies
           write(iOIO,30) lsp(i), wpar(i), nsp(i,1), nsp(i,2)
        enddo
        write(iOIO,*) ' '
      endif
#endif
      
 30   format (i10,2x,g12.5,2x,2(i8,2x))
CEVERINO112205: It was above #ifdef STARFORM
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
c      
      close (iOIO)
c
      return
      end
c
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
c     
c     ------------------------------------------
      subroutine Read_Particle_Time ( FileName )
c     ------------------------------------------
c
c     purpose: when the code starts or restarts
c              this routine takes care of array pt
c 
      include 'a_tree.h'
      include 'a_control.h'
      character*256 FileName 
c
      open ( 60 ,file = TRIM(FileName),
     &       form = 'unformatted',status = 'old' )
      read (60) (pdt(ic1),ic1=1,np)
      close( 60 )
c
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
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'
      character*256 FileName 
      real*8 tdum, adum 
      real*8 ws_old, ws_oldi
      common / s_old1 / ws_old, ws_oldi
c
      open ( 60 ,file = TRIM(FileName),
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
