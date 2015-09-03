c
c     Analysis of cosmological simulations
c
c     version which analyze multiples snatshoots.
c
c     by Daniel Ceverino, Jerusalem, 2010
c
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      include 'analysis_Rec.h'
      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' /
c
      call ReadControl ()
c
      call ReadSchedule ()
c
      if ( iflag0d.eq.1) then
         write (*,*)  'Setting units...'
         call Create_File_Names(0)
         call Read_Gas_Binary ( dt , fname_gas  )
         call units()           ! Set physical units
         write (*,*)  'Setting global variables...'
         rho0C =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) ! (comoving) n_H in cm^-3
         gamma1 = gamma -1.00
         do Level=MinLevel, MaxLevel
            CellSize(Level)=1.0 / 2**Level
         enddo
         lprefix=index(prefix_name, ' ') -1
      endif

c
c     Set output Files for evolution analysis:
c

c
      DO i=1,n_save ! loop over files
         aexpn=asave(i) 
         ic = int(aexpn*1000)
         i1 = ic / 1000 
         i2 = (ic - (ic/1000)*1000)/100 
         i3 = (ic - (ic/100)*100) / 10 
         i4 = (ic - (ic/10)*10)
         sufix_name='a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
         lsufix = index(sufix_name    , ' ') - 1
         write (*,*)  '--------------------------'
         write (*,*)  'CURRENT SNAPSHOT: ', aexpn
         write (*,*)  '--------------------------'
c
         call Init_Arrays           ()
c
c     READING FILES:
         call Create_File_Names(2)
         if ( iflag0a.eq.1) then
            write (*,11)    fname_part1
            write (*,11)    fname_part2         
            call Read_Particles_Binary(fname_part1,fname_part2)
         endif
         if ( iflag0b.eq.1) then
            write (*,11)    fname_part4
            call   Read_Stellar_Data ( fname_part4)
         endif
         if ( iflag0c.eq.1) then
            write (*,11)   fname_gas
            call Read_Gas_Binary ( dt , fname_gas  )
         endif
c-------------------------------------------------------------------- 
c     PRE-ANALYSIS OF EACH FILE:
c-------------------------------------------------------------------- 
         rminC  =rminC0/ ( r0 * 1000.0 / hubble * aexpn) ! in code units
         rmaxC  =rmaxC0/ ( r0 * 1000.0 / hubble * aexpn)
         Rs =Rs0 / ( r0 * 1000.0 / hubble * aexpn) ! in code units
      if (iflag0e.eq.1) then     !! For particles
         write (*,*)  'Center position and velocity of particles:'
         write (*,*) '-------------------------------------------------'
         xc=(boxh/2.)/r0 ! Initial center in code units.
         yc=(boxh/2.)/r0
         zc=(boxh/2.)/r0
         if ( ioptCenter.eq.1) then
            call Find_Halo_Center ( xc, yc, zc, 3 ) ! potential
            write (*,*) 'Position of minimum potential in  proper Kpc:'
         elseif ( ioptCenter.eq.2) then
            write (*,*)'Initial center suplied by user'
            write (*,*)'from Halo finder data'
            xc = (xuser  / r0 ) +1. ! From Halo Finder
            yc = (yuser  / r0 ) +1.
            zc = (zuser  / r0 ) +1.
         elseif ( ioptCenter.eq.3) then
            write (*,*)'Initial center suplied by user'
            write (*,*)'from IFRIT (original), '
            write (*,*)  'the code takes care about the flip of axis'
            xc=(zuser +1.)/2. * boxh / r0  +1 ! Be careful with the flip of axis in ifrit!!
            yc=(yuser +1.)/2. * boxh / r0  +1
            zc=(xuser +1.)/2. * boxh / r0  +1
         endif
         write (*,*)  (xc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)  (yc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)  (zc-1.)*1000.0*aexpn*r0/hubble
c----ST: removed starcenter for DM only:
         call STARcenter (rminC, rmaxC, Xc , Yc , Zc, Vxc, Vyc, Vzc)
c----------------
         write (*,*) 'Center in proper Kpc:'
         write (*,*)   (xc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)   (yc-1.)*1000.0*aexpn*r0/hubble 
         write (*,*)   (zc-1.)*1000.0*aexpn*r0/hubble
         write (*,*)  'Velocity of the center in Km/s'
         write (*,*)  v0*Vxc / aexpn
         write (*,*)  v0*Vyc / aexpn
         write (*,*)  v0*Vzc / aexpn
      endif
!     Find total angular momentum:
      if (iflag0f.eq.1) then
         r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
         v0P = v0 / aexpn       ! in km/s (proper units)
         if (isys.eq.3) then
            Tgas=1.e4
          write (*,*) 'Computing angular momentum for cold (T<10^4) gas'
         else
            Tgas=1.e9
         endif
         call angularM(Vxc,Vyc,Vzc,xc,yc,zc,isys,Rs,Tgas,aLx,aLy,aLz)
         write (*,*)  'Angular momentum in Msun Kpc Km/s (proper):'
         write (*,14)  aLx*aM0*r0P*v0P
         write (*,14)  aLy*aM0*r0P*v0P
         write (*,14)  aLz*aM0*r0P*v0P
 14      format(g12.3)
         aL = sqrt ( aLx**2 + aLy**2 + aLz**2)
         write (*,*)  'AM axis'
         write (*,*)  aLx/aL
         write (*,*)  aLy/aL
         write (*,*)  aLz/aL
         call NewAxes(aLx, aLy, aLz, AxX,AxY, AxZ)
         write (*,*) ' Vector 1= ', AXX(1), AxX(2), AxX(3)
         write (*,*)  'Vector 2= ', AXY(1), AXY(2), AXY(3)
         write (*,*)  'Vector 3= ', AxZ(1), AxZ(2), AxZ(3)
      endif
c-------------------------------------------------------------------- 
c     ANALYSIS OF EACH FILE:
c-------------------------------------------------------------------- 
      rmin=rmin0/ ( r0 * 1000.0 / hubble * aexpn)
      rmax=rmax0/ ( r0 * 1000.0 / hubble * aexpn)
      ZDmax=ZDmax0/( r0 * 1000.0 / hubble * aexpn)
      if (iflag1.eq.1) then  
         write (*,*) " writing an ASCII file with"
         write (*,340) tAge
         write (*,*)  'in spherical coordinates'
         write (*,*)  'inside 5 Kpc'
         rmaxASCII=5.0/ ( r0 * 1000.0 / hubble * aexpn)
         filename='ASCII_St_Test' // sufix_name(1:lsufix)
            call Write_Stars_ASCII(filename, Vxc, Vyc, Vzc,xc,yc,zc, 
     &             rmaxASCII,tAge,ZDmax,AxX,AxY,AxZ,1)
 340  format ( 'stars younger than ', E6.1, ' yrs')
      endif
      if (iflag2.eq.1) then ! Spherical profiles of density, temperature,..
         filename='./profileSpGe_Test' // sufix_name(1:lsufix)
         call ProfilesSp (xc,yc,zc,T_gas,rmin, rmax,nrbin,filename,
     &                    Rvir,aMvir)
         write (*,*)  'Rvir= ', Rvir ,' kpc'
         write (*,*)  'Mvir=' , aMvir, ' Msun'
         rmax= 2 * Rvir / ( r0 * 1000.0 / hubble * aexpn)
         rmax2= 4 * Rvir 
         if (rmax2.lt.100.) then
            write (middle_name, 341) rmax2
            middle_name='0'//middle_name
 341     format (F3.0)
         else
            write (middle_name, 342) rmax2
 342        format (F4.0)
         endif
         
ccccccccccc Kenza for Jacob
         call Open_ASCII_File(100, 'summary.txt', append)
c         write (100,*) 'aexpn	redshift  Mvir        Rvir       a11  ',
c     &   'a12           a13         a21         a22         a23       ',
c     &   'a31         a32           a33         vpec_x  vpec_y  vpec_z'

         write (100,61) 'aexpn', 'redshift', 'Mvir', 'Rvir', 'a11',
     &   'a12', 'a13', 'a21', 'a22', 'a23', 'a31', 'a32', 'a33', 
     &   'vpec_x', 'vpec_y', 'vpec_z'
 61      format (A5, 5x, A8, 5x, A4, 9x, A4, 9x,  A3, 11x, 
     &        A3, 11x,  A3, 11x,  A3, 11x,  A3, 11x,  A3, 11x,  A3, 11x, 
     &        A3, 11x,  A3, 9x,  A6, 9x,  A6, 9x,  A6)
         write (100,62) '(M_solar)', '(kpc)', '(km/s)', 
     &        '(km/s)', '(km/s)'
 62      format(23x, A9, 4x, A5, 132x, A6, 9x, A6, 9x, A6)
         write (100,65)  aexpn, (1./aexpn - 1.), aMvir, Rvir, 
     &        AXX(1), AxX(2), AxX(3), AXY(1), AXY(2), 
     &        AXY(3), AxZ(1), AxZ(2), AxZ(3), 
     &        v0*Vxc / aexpn, v0*Vyc / aexpn,  v0*Vzc / aexpn
 65      format(F10.6, 5x, F8.5, 3x, E12.4, 3x, F8.4, 4x, 
     &        F11.7, 2x, F11.7, 3x, F11.7, 3x, F11.7, 3x, F11.7, 3x, 
     &        F11.7, 3x, F11.7, 3x, F11.7, 3x, F11.7, 3x, 
     &        F11.5, 5x, F11.5, 5x, F11.5)

         close (100)

c========STG:
         write (middle_name, '(a)') '' 
c===========
         lmiddle= index( middle_name, ' ' ) -1
      endif
      write (*,*)  'Writing binary files..'
      if (iflag3.eq.1) then
         filename=  
     &   prefix_name(1:lprefix)// '_D'// middle_name(1:lmiddle) //
     &   sufix_name(1:lsufix) 
         call Write_DM_binary(filename,Vxc,Vyc,Vzc,xc,yc,zc,rmax,
     &              AxX,AxY,AxZ,0)
      endif
      if (iflag4.eq.1) then
         filename=  
     &   prefix_name(1:lprefix)// '_S'// middle_name(1:lmiddle) // 
     &   sufix_name(1:lsufix) 
         write (*,340) tAge
         call Write_Stars_binary(filename,Vxc, Vyc, Vzc,xc,yc,zc, 
     &                 rmax,tAge,ZDmax,AxX,AxY,AxZ,0)
      endif
      if (iflag5.eq.1) then
         filename=  
     &   prefix_name(1:lprefix)// '_Si'// middle_name(1:lmiddle) //
     &   sufix_name(1:lsufix) 
         write (*,340) tAge
         call Write_Stars_binary0(filename,Vxc, Vyc, Vzc,xc,yc,zc, 
     &                 rmax,tAge,ZDmax,AxX,AxY,AxZ,0)
      endif
      if (iflag6.eq.1) then
         filename=  
     &   prefix_name(1:lprefix)// '_SZ'// middle_name(1:lmiddle) //
     &   sufix_name(1:lsufix) 
         write (*,340) tAge
         call Write_Stars_binaryZ(filename,Vxc, Vyc, Vzc,xc,yc,zc, 
     &                 rmax,tAge,ZDmax,AxX,AxY,AxZ,0)
      endif
      if (iflag7.eq.1) then
         Ext_box = 4 * Rvir
         filename=  
     &   prefix_name(1:lprefix)// '_G'// middle_name(1:lmiddle) //
     &   sufix_name(1:lsufix)
         call Extract_BOX3(filename,xc,yc,zc,Vxc,Vyc,Vzc,Ext_box)
      endif
      if (iflag8.eq.1) then
         Ext_box = 4 * Rvir
         filename=  
     &   prefix_name(1:lprefix)// '_GZ'// middle_name(1:lmiddle) //
     &   sufix_name(1:lsufix)
         call Extract_BOX6(filename,xc,yc,zc,Vxc,Vyc,Vzc,Ext_box)
      endif

      ENDDO
 11   format(  'Reading  file: ', A30)
      end
c-------------------------------------------------------------------- 
      SUBROUTINE ReadSchedule ()
c-------------------------------------------------------------------- 
c     Read the file schedule_R.dat
c
      include 'analysis_Rec.h'
      character*60 title
c
      OPEN (63, file='schedule_R.dat')
C-----------------Schedule of analysis (flags= 0 Off , 1 On)-------- 
      read (63,10) title
 10   format (A60)
      write (*,10) title

      read (63,*) iflag0a
      if ( iflag0a.eq.1) write (*,*) 'Read N-Body file'
      read (63,*) iflag0b
      if ( iflag0b.eq.1) write (*,*) 'Read stellar file'
      read (63,*) iflag0c
      if ( iflag0c.eq.1) write (*,*) 'Read HYDRO file'
      read (63,*) iflag0d
      if ( iflag0d.eq.1) write (*,*) 'Set units'
      read (63,*) iflag0e
      if ( iflag0e.eq.1) write (*,*)  'Find center using particles'
      read (63,*) iflag0f
      if (iflag0f.eq.1) write (*,*) 'Compute the angular momentum & 
     & new system of coordinates'
      read (63,*) iflag1
      if (iflag1.eq.1)write (*,*)'ASCII of first specie of DM or stars'
      read (63,*) iflag2
      if (iflag2.eq.1) write (*,*) 'spherical profiles'
      read (63,*) iflag3
      if (iflag3.eq.1) 
     & write (*,*) 'Binary file with DM particles inside 4Rvir-box'
      read (63,*) iflag4
      if (iflag4.eq.1) 
     & write (*,*) 'Binary file with stellar particles inside 4Rvir box'
      read (63,*) iflag5
      if (iflag5.eq.1) 
     & write (*,*) 
     &'Binary file with stellar particles inside 4Rvir (initial mass)'
      read (63,*) iflag6
      if (iflag6.eq.1) 
     & write (*,*) 
     &'Binary file with stellar particles inside 4Rvir (metallicities)'
      read (63,*) iflag7
      if(iflag7.eq.1)
     & write(*,*)'Binary file with gas cells inside 4Rvir box'
      read (63,*) iflag8
      if(iflag8.eq.1)
     & write(*,*)'Binary file with gas cells inside 4Rvir (metals)'

c-------------------------------------------------------------------- 
c    Parameters for PRE-analysis
C-------------------------------------------------------------------- 
      read (63,10) title
      read (63,*) gamma
!     Find center of the main Halo:    
      read (63,*) rminC0  ! in proper Kpc
      read (63,*) rmaxC0 
      read (63,*) ioptCenter     ! =1/2/3--> potential/HF/IFRIT
      read (63,*) xuser          ! if ioptCenter=2/3
      read (63,*) yuser
      read (63,*) zuser
      read (63,*) isys           ! For angular momentum calculations:  isys=1 --> DM particles,  isys=2 --> stars, isys=3 --> cold gas
      read (63,*) Rs0             ! in proper Kpc
c-------------------------------------------------------------------- 
c    Parameters for Analysis
C-------------------------------------------------------------------- 
      read (63,10) title
      read (63,*) rmin0           ! min. radius of the profile
      read (63,*) rmax0           ! max. radius of the profile
      read (63,*) nrbin          ! number of bins
      read (63,20) prefix_name  
      read (63,*) icomponent ! 0 / 1 / 2 / 3 --> all / gas /stars / dm
      read (63,*) T_gas
      read (63,*) tAge
      read (63,*) Zdmax0
 20   format (A256)


      end
c-------------------------------------------------------------------- 
      SUBROUTINE  Create_File_Names(iMode)
c-------------------------------------------------------------------- 
c
c     iMode defines how the output file is named
c           = 0 , filename is just the jobname
c           = 1 , filename is jobname + global step number
c           = 2 , filename is jobname + expansion parameter       
c
#     include "a_def.h"
      include 'a_setup.h'
      include 'a_control.h'
      include 'analysis_Rec.h'

      character digits(0:9)
      data digits / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' / 
      character*5 chstep
      real*8 dt

      ljob  = index(jobname1, ' ') - 1
      lpath = index(path    , ' ') - 1

      if ( iMode .eq. 0 ) then 
        fname_gas  = 
     &  path(1:lpath)//'/'//jobname1(1:ljob)//'.d '
        fname_part1 = 
     &  path(1:lpath)//'/PMcrd.DAT '
        fname_part2 = 
     &  path(1:lpath)//'/PMcrs0.DAT '
        fname_part3 = 
     &  path(1:lpath)//'/pt.dat '
        fname_part4 = 
     &  path(1:lpath)//'/stars.dat '
      endif
      if ( iMode .eq. 1 ) then
        is = istep 
        idummy = is / 10000
        is = is - idummy * 10000
        chstep = digits(idummy)
        ich = 1
        do i = 3 , 0 , -1
          idummy = is / 10**i        
          chstep = chstep(1:ich) // digits(idummy)
          ich = ich + 1
          is = is - idummy*10**i
        enddo
        fname_gas  = 
     &  path(1:lpath)//'/'//jobname1(1:ljob)//'_'//chstep//'.d '
        fname_part1 = 
     &  path(1:lpath)//'/PMcrd_'//chstep//'.DAT '
        fname_part2 = 
     &  path(1:lpath)//'/PMcrs0_'//chstep//'.DAT '
        fname_part3 = 
     &  path(1:lpath)//'/pt_'//chstep//'.dat '
        fname_part4 = 
     &  path(1:lpath)//'/stars_'//chstep//'.dat '
      endif
      if ( iMode .eq. 2 ) then 
          ic = int(aexpn*1000)
          i1 = ic / 1000 
          i2 = (ic - (ic/1000)*1000)/100 
          i3 = (ic - (ic/100)*100) / 10 
          i4 = (ic - (ic/10)*10)

          fname_gas = path(1:lpath)//'/'//jobname1(1:ljob)//'_a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.d '
          fname_part1 = path(1:lpath)//'/PMcrda'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.DAT '
          fname_part2 = path(1:lpath)//'/PMcrs0a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.DAT '
          fname_part3 = path(1:lpath)//'/pta'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
          fname_part4 = path(1:lpath)//'/stars_a'//
     &            digits(i1)//'.'//digits(i2)//digits(i3)//digits(i4)
     &            //'.dat '
       endif
       end
