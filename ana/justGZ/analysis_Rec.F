c
c     Analysis of cosmological simulations
c     version which analyzes multiple snapshots
c     by Daniel Ceverino, Jerusalem, 2010
c
c     edits by Kenza Arraki, 2014
c
      include 'a_tree.h'
      include 'a_control.h'
c
      character*5 aexpnstr
      character*256 outpath
      character*256 fname_gas, fname_part1, fname_part2,
     &               fname_part3, fname_part4, fname_read,
     &               filename
      real*8 xuser, yuser, zuser, Rs0, Rvir_i0, vmax, Mvir, Rvir
      integer Npb
      real*8 xc, yc, zc, Rs, rmin, rmax, rmax2
      real*8 Vxc, Vyc, Vzc
      real*8 Mtot, Mdm, Mgas, Mst, Tstars
      real*8 AxX(3), AxY(3), AxZ(3) 
      real*8 aLx, aLy, aLz, aL
      integer i, j
      integer, parameter :: nrbin = 100     !# bins for analysis
      real*8, parameter :: Tgas_prof = 1.0e9 !max gas temp for profiles
      real*8, parameter :: Tgas_angm = 1.0e4 !max gas temp for angular momentum
      real*8, parameter :: tAge  = 1.0e11  !max stellar age for all outputs
      real*8, parameter :: rmin0 = 0.01    !min radius, kpc (~force softening)
      outpath = 'output'
      WRITE (iofile, '(a,a)') TRIM(outpath),'/ioread.txt'
c
c     
      write (*,'(/1x, ''****************************************'',
     &           /1x, ''* ANALYSIS of hydroART ouputs  -  2014 *'',
     &           /1x, ''****************************************'')')
c
      call ReadControl ()
c
c
cc    LOOP OVER ALL PROVIDED FILES:
      DO i=1,n_save
         aexpn=asave(i) 
         if(aexpn.lt.1)then
            WRITE (aexpnstr,'(i1.1,a,i3.3)') INT(aexpn),'.',
     &        MIN(INT((aexpn)*1000.+0.5),999)
         else
            WRITE (aexpnstr,'(i1.1,a,i3.3)') INT(aexpn),'.',
     &        MIN(INT((aexpn-1.)*1000.+0.5),999)
         endif         
         write (*,*)  '--------------------------'
         write (*,*)  '--------------------------'
         write (*,*)  'CURRENT SNAPSHOT: ', aexpn
         call Init_Arrays ()
c     
cc       READING FILES:
         WRITE(fname_gas,'(a,a,a,a,a,a)') TRIM(path),'/',
     &        TRIM(jobname1),'_a',aexpnstr,'.d'
         WRITE(fname_part1,'(a,a,a,a)') TRIM(path),'/PMcrda',
     &        aexpnstr,'.DAT'
         WRITE(fname_part2,'(a,a,a,a)') TRIM(path),'/PMcrs0a',
     &        aexpnstr,'.DAT'
         WRITE(fname_part3,'(a,a,a,a)') TRIM(path),'/pta',
     &        aexpnstr,'.dat'
         WRITE(fname_part4,'(a,a,a,a)') TRIM(path),'/stars_a',
     &        aexpnstr,'.dat'
c
         write (*,*)  '--------------------------'
         write(*,'(a)') 'Reading  files: '
         write(*,'(a)') TRIM(fname_gas)
         call Read_Gas_Binary (dt, fname_gas)
         write(*,'(a)') TRIM(fname_part1), TRIM(fname_part2)
         call Read_Particles_Binary(fname_part1,fname_part2)
         write(*,'(a)') TRIM(fname_part4)
         call Read_Stellar_Data (fname_part4)
c
cc       SETTING UNITS:
         write (*,*)  'Setting units.'
         call SetUnits()        ! Set physical units
         rho0C =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) ! comoving n_H in cm^-3
         do Level=MinLevel, MaxLevel
            CellSize(Level) = 1.0/2.0**Level
         enddo
         r0P = r0*1000.0*aexpn/hubble ! in physical kpc
         v0P = v0/aexpn         ! in proper km/s
c
cc       DO ANALYSIS
         write (*,*) '--------------------------'
         write (*,*) '--------------------------'
         write (*,*) '       Main Halo'
         !create input file name and read the first line
         WRITE (fname_read, '(a,a,a,a)') TRIM(outpath),'/input_',
     &        aexpnstr,'.txt'
         open (63, file=TRIM(fname_read))
         call ReadInput(fname_read,xuser,yuser,zuser,Rs0,Rvir_i0,vmax,
     &        Npb)
         close (63)
c
cc       PRE-ANALYSIS OF EACH FILE:
         xc = (xuser / r0) + 1. ! center from Halo Finder
         yc = (yuser / r0) + 1. ! /h comoving Mpc to code units
         zc = (zuser / r0) + 1.
         Rs = Rs0/r0P           ! physical kpc to code units
         Rvir = Rvir_i0         ! virial radius in physical kpc
         rmin = rmin0/r0P            
         rmax = Rvir_i0/r0P     ! virial radius in code units
c                               ! vmax is in km/s
         write (*,*) '--------------------------'
         write (*,*) 'Center in /h comoving Mpc, proper kpc:'
         write (*,*) (xc-1.)*r0, (xc-1.)*r0P
         write (*,*) (yc-1.)*r0, (yc-1.)*r0P
         write (*,*) (zc-1.)*r0, (zc-1.)*r0P
c     
cc       FIND HALO CENTER VELOCITY
         write (*,*) '--------------------------'
         write (*,*) 'Velocity of the center in km/s:'
         call Find_Halo_Velocity (Rs, xc, yc, zc, Vxc, Vyc, Vzc)
         write (*,*) Vxc*v0P
         write (*,*) Vyc*v0P
         write (*,*) Vzc*v0P
c
cc       ANALYSIS
         write (*,*) 'Writing halo profile'
         WRITE (filename, '(6(a))') TRIM(outpath),'/allprofile_a',
     &        aexpnstr,'_1.txt'
         Mtot = 0.
         Mdm = 0.
         Mgas = 0.
         Mst = 0.
         Tstars = 0.
         call AllProfile (filename,xc,yc,zc,rmin,rmax,nrbin,
     &        Mtot,Mdm,Mgas,Mst,Tstars)
         Mvir = Mtot
         write (*,'(a,f12.4,a)') ' Rvir = ', Rvir,' kpc'
         write (*,'(a,e12.4,a)') ' Mvir = ', Mvir,' Msun'
c
cc       TOTAL ANGULAR MOMENTUM using stars and DM (j=1,2)
         do j=1,3
            write (*,*) '----------------------------------------'
            call angularM(xc,yc,zc,Vxc,Vyc,Vzc,j,Rs,
     &           Tgas_angm,aLx,aLy,aLz)
            write (*,*) 
     &           'Angular momentum in proper Msun kpc km/s:'
            write (*,'(g12.3)') aLx*aM0*r0P*v0P
            write (*,'(g12.3)') aLy*aM0*r0P*v0P
            write (*,'(g12.3)') aLz*aM0*r0P*v0P
            aL = sqrt(aLx**2 + aLy**2 + aLz**2)
            write (*,*)  'Angular Momentum axes:'
            write (*,'(f12.4)') aLx/aL
            write (*,'(f12.4)') aLy/aL
            write (*,'(f12.4)') aLz/aL
cc          FIND TOTAL ANGULAR MOMENTUM of cold gas (j=3)
            if (j .eq. 3) then
               ! Find a new system of coordinates
               call NewAxes(aLx,aLy,aLz,AxX,AxY,AxZ)
               write (*,*) 'New axes:'
               write (*,'(a,3(g12.4,1x))') ' Vector 1 = ',
     &              AxX(1),AxX(2),AxX(3)
               write (*,'(a,3(g12.4,1x))') ' Vector 2 = ',
     &              AxY(1),AxY(2),AxY(3)
               write (*,'(a,3(g12.4,1x))') ' Vector 3 = ',
     &              AxZ(1),AxZ(2),AxZ(3)
            endif
         enddo
c     
cc       write out specific information for QSO studies
         WRITE (filename, '(5(a))') TRIM(outpath),'/rotmat_',
     &        'a',aexpnstr,'.txt'
         open (64, file=TRIM(filename))
         write(64, '(3(a))') 'aexpn redshift Mvir Rvir ',
     &        'a11 a12 a13 a21 a22 a23 a31 a32 a33 ',
     &        'vpec x  vpec y  vpec z'
         write (64, '(a,15(g15.6))')
     &        aexpnstr, (1./aexpn)-1., Mvir, Rvir, AxX, AxY, AxZ,
     &        Vxc*v0P, Vyc*v0P, Vzc*v0P
         close (64)
c
         write (*,*) '-----------------------------------------'
         write (*,*)  'Writing zoom files..'
         ! Write binary/ASCII files zoomed in on the selected halo
         ! centered on a halo and out to a given radius, rmax
         ! The output box has rmax = L/2, so L=2*rmax
c        DM FILE
         WRITE(filename,'(6(a))') TRIM(outpath),'/',TRIM(jobname1),
     &        '_DMa',aexpnstr,'.dat'
         rmax2 = 3.*Rvir/r0P
         call Write_DM_binary_zoom(filename,xc,yc,zc,Vxc,Vyc,Vzc,
     &        rmax2)
c        STELLAR FILE
         WRITE(filename,'(6(a))') TRIM(outpath),'/',TRIM(jobname1),
     &        '_SZa',aexpnstr,'.dat'
         rmax2 = 3.*Rvir/r0P
         call Write_Stars_binary_zoom(filename,xc,yc,zc,Vxc,Vyc,Vzc,
     &        rmax2)
c        GAS FILE (as ASCII)
         WRITE(filename,'(6(a))') TRIM(outpath),'/',TRIM(jobname1),
     &        '_GZa',aexpnstr,'.txt'
         rmax2 = 3.*Rvir/r0P
         call Write_Gas_ASCII_zoom(filename,xc,yc,zc,Vxc,Vyc,Vzc,
     &        rmax2)

ccccc TEST ALTERNATE GAS FILE WRITE
c         WRITE(filename,'(6(a))') TRIM(outpath),'/',TRIM(jobname1),
c     &        '_GZa',aexpnstr,'_2.dat'
c         call Extract_BOX6(filename,xc,yc,zc,Vxc,Vyc,Vzc,4.*Rvir)
c
      ENDDO
      end
c      
c-------------------------------------------------------------------- 
      SUBROUTINE ReadInput(filename,xuser,yuser,zuser,Rs0,Rvir_i0,
     &     vmax,Npb)
cc     Read the input file
c      x center, y center, z center (all in comoving /h Mpc), scale radius (proper kpc), approximate virial radius (proper kpc), v_max (km/s), # bound particles
c
      real*8 xuser, yuser, zuser, Rs0, Rvir_i0, vmax
      integer Npb
      character*256 filename
      xuser=0.
      yuser=0.
      zuser=0.
      Rs0=0.
      Rvir_i0=0.
      vmax=0.
      Npb=0
c
      read (63,*) xuser, yuser, zuser, Rs0, Rvir_i0, vmax, Npb
c
      write (*,*) 'Parameters from input.txt'
      write (*,*) xuser  
      write (*,*) yuser
      write (*,*) zuser
      write (*,*) Rs0    
      write (*,*) Rvir_i0  
      write (*,*) vmax
      write (*,*) Npb
      end

