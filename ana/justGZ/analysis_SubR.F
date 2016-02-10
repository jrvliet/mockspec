
c   Subroutines for analysis of cosmological simulations ( recursive version)
c     by Daniel Ceverino, Jerusalem (2010)
c
c     edited by Kenza Arraki, 2014
c
c
c-------------------------------------------------------------------- 
      SUBROUTINE Write_ASCII_Particles(filename)
c     Write out all particles in ASCII format
c     Designed for Rockstar-Galaxies to readin with correct format
c     Writes:   X  Y  Z  VX  VY  VZ  MASS  ENERGY  ID
c     X Y and Z in comoving Mpc/h (cMpc/h)
c     VX VY and VZ in km/s
c     Mass in Msun/h   ENERGY set to 0   ID set to integer particle #
c
      include 'a_tree.h'
      include 'a_control.h'
c
      character*256 filename
      integer io, ip
      io = 40
c
      open(io, file=TRIM(filename), form='formatted')
c
      write(io,'(a)') '#Write out all particles in ASCII format'
      write(io,'(a)') '#X  Y  Z  VX  VY  VZ  MASS  ENERGY  ID'
      write(io,'(a)') '#comoving /h Mpc; km/s; /h Msun; energy=0'
      write(io,'(a,a)') '#file = ',TRIM(jobname1)
      write(io,'(a,f12.6)') '#a = ',aexpn
      write(io,'(a,f12.6)') '#h = ',hubble
      write(io,'(a,f12.6)') '#Ol = ',Oml0
      write(io,'(a,f12.6)') '#Om = ',Om0
      write(io,'(a)') '# '
      do ip = 1, nsp(nspecies,2)
         write(io,'(8(e16.6),i)') (x(ip)-1.)*r0, (y(ip)-1.)*r0, 
     &      (z(ip)-1.)*r0, vx(ip)*v0/aexpn, vy(ip)*v0/aexpn, 
     &      vz(ip)*v0/aexpn, pw(ip)*aM0*hubble, 0.0, ip
      enddo
      close(io)
      return
      END
c-------------------------------------------------------------------- 
      SUBROUTINE Check_Masses ()
c-------------------------------------------------------------------- 
c     check the mass of all the particles and cells
c
c     input in code units
      include 'a_tree.h'
      include 'a_control.h'
c
      integer nn, icell0, icell, Level, nLevel, i1, i2, idc
      real*8 cellvolume, rhototal, mstartot, mdmtot
c
      write(*,*) nspecies
      write(*,'(e16.7)') pw(1)*aM0
      write(*,'(e16.7)') pw(nsp(1,1))*aM0
      write(*,'(e16.7)') pw(nsp(2,1))*aM0
      write(*,'(e16.7)') pw(nsp(3,1))*aM0
      write(*,'(e16.7)') pw(nsp(4,1))*aM0
      write(*,'(e16.7)') pw(nsp(5,1))*aM0
c
cc     sum up the mass in all particles
      write(*,'(a,g16.7,a)') 'all particles', SUM(pw)*aM0, 'Msun'
c      sum up each species
      mdmtot = 0.
      do ip=1,nsp(nspecies-1,2)
         mparttot = mparttot+pw(ip)
      enddo
      write(*,'(a,g16.7,a)') 'darkmatter', mdmtot*aM0, 'Msun'
      mdmtot = 0.
      do ip=nsp(1,1),nsp(1,2)
         mdmtot = mdmtot+pw(ip)
      enddo
      write(*,'(a,g16.7,a)') 'specie0', mdmtot*aM0, 'Msun'
      mdmtot = 0.
      do ip=nsp(2,1),nsp(2,2)
         mdmtot = mdmtot+pw(ip)
      enddo
      write(*,'(a,g16.7,a)') 'specie1', mdmtot*aM0, 'Msun'
      mdmtot = 0.
      do ip=nsp(3,1),nsp(3,2)
         mdmtot = mdmtot+pw(ip)
      enddo
      write(*,'(a,g16.7,a)') 'specie2', mdmtot*aM0, 'Msun'
      mdmtot = 0.
      do ip=nsp(4,1),nsp(4,2)
         mdmtot = mdmtot+pw(ip)
      enddo
      write(*,'(a,g16.7,a)') 'specie3', mdmtot*aM0, 'Msun'
c
c
      mstartot = 0.
      do ip=nsp(nspecies,1),nsp(nspecies,2)
         mstartot = mstartot+pw(ip)
      enddo
      write(*,'(a,g16.7,a)') 'specie4', mstartot*aM0, 'Msun'
      write(*,'(a,g16.7,a)') 'stars', SUM(pw0)*aM0, 'Msun'
c
cc
      nn = 0
      rhototal = 0.
      call Get_MaxLevelNow ()
      do Level= MinLevel, MaxLevelNow
         if (Level .eq. MinLevel) then
            cellvolume = CellSize(Level)**3
            do icell0=1, ncell0
               rhototal = rhototal + hvar(1,icell0) * cellvolume
               nn = nn+1
            enddo
c            write(*,'(a,g12.4,a)') 'root mesh rho',rhototal*rho0,'Msun'
         else
            cellvolume = CellSize(Level)**3
            call Select_Cells ( Level , nLevel )
            do i1 = 1 , nLevel
               icell = iSelect(i1)
               do i2 = 0, 7
                  idc = icell + i2
                  rhototal = rhototal + hvar(1,idc) * cellvolume
                  nn = nn+1
               enddo
            enddo
         endif
         write(*,*) Level, nn
      enddo
      write(*,'(a,g12.4,a)') 'total rho', rhototal*aM0, 'Msun'
      return
      end
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      subroutine Find_Halo_Center_close (xc, yc, zc, rmax)
c     Finding a new center of the halo
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 xc, yc, zc, rmax
c
      rhogmax = -1.e6           ! max. gas density 
      rhotmax = -1.e6           ! max. total density
      phimin  = 1.e6            ! minimum of potential
      rhodmmax = -1.e6 
      irhogmax = -1
      irhotmax = -1 
      iphimin  = -1 
      irhodmmax = -1
c     
      call Get_MaxLevelNow ()
      DO Level = MinLevel+1, MaxLevelNow  
         call Select_Cells ( Level , nLevel ) 
         do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
               ic = icell + ic2
               if ( iOctCh(ic) .eq. nil ) then 
                  call Ps ( ic , Posx , Posy , Posz )
c                  call PeriodicBoundary( xc, yc, zc, Posx, Posy, Posz)
                  rr = sqrt((Posx-xc)**2+(Posy-yc)**2+(Posz-zc)**2)
                  if ( rr .lt. rmax ) then 
                     dcell1 = var(1,ic)*2.0**(3.*Level) + 1.0 ! tot. density
                     dcell2 = var(2,ic) ! potential 
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
                  endif
               endif
            enddo
         enddo
      ENDDO   
c
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
c     
      return
      end
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE Find_Halo_Velocity (Rs, xc, yc, zc, Vxc, Vyc, Vzc)
c-------------------------------------------------------------------- 
c     find the velocity of a halo given its center position
c
c     input in code units
      include 'a_tree.h'
      include 'a_control.h'
c
      integer nn, ip
      real*8 Rs, Rselect, R2, rp2
      real*8 xc, yc, zc, Vxc, Vyc, Vzc 
      real*8 v8x, v8y, v8z
      integer, parameter :: maxnn = 5000, minnn = 20

      Rselect = Rs
      nn = 0
c
cc     find the average velocity of the center particles
cc     keep the results if the # particles is between min and max
      do while ((nn .gt. maxnn) .or. (nn .lt. minnn))
         R2  = Rselect**2
         nn  = 0
         v8x = 0.
         v8y = 0.
         v8z = 0.
         do ip = 1, nsp(nspec,2) !loop over all particles
            rp2 = (x(ip) - xc)**2 + (y(ip) - yc)**2 + (z(ip) - zc)**2
            if (rp2 .lt. R2) then !if particle closer than R2
               nn = nn + 1
               v8x = v8x + vx(ip)
               v8y = v8y + vy(ip)
               v8z = v8z + vz(ip)
            endif
         enddo
         Vxc = v8x / nn
         Vyc = v8y / nn
         Vzc = v8z / nn
         if(nn .lt. 20) then
            Rselect = Rselect * 1.1
         else
            Rselect = Rselect / 1.1
         endif
      enddo     
      write (*,'(a,i4,a,f6.3,a)') 'Calculated using: ',
     &  nn, ' particles within ',
     & (Rselect)*1000.0*aexpn*r0/hubble, ' proper kpc'
      return
      end
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE Write_some_ASCII_Particles(filename,xc,yc,zc,rmax)
c     write out particles that are close to the main halo
c
c     input in code units
      include 'a_tree.h'
      include 'a_control.h'
#     include 'a_def.h'
c
      character*256 filename
      real*8 xc, yc, zc, rmax, rpart, rmax2
      integer ip
      integer io
      io = 41
c
      open(io, file=TRIM(filename), form='formatted')
      write(io,'(a)') '#Write out all particles in ASCII format'
      write(io,'(a)') '#X  Y  Z  VX  VY  VZ  MASS  ENERGY  ID'
      write(io,'(a)') '#comoving /h Mpc; km/s; /h Msun; energy=0'
      write(io,'(a,a)') '# file = ',TRIM(jobname1)
      write(io,'(a,f12.6)') '# a = ',aexpn
      write(io,'(3(a,f12.6))')
     &     '# h = ',hubble,' Ol = ',Oml0,' Om = ',Om0
      write(io,'(a)') '# '
c
      rmax2 = rmax**2
      do ip = nsp(1,1), nsp(nspecies,2) !for all particles
         rpart = (x(ip)-xc)**2 + (y(ip)-yc)**2 + (z(ip)-zc)**2
         if (rpart .le. rmax2) then
            write(io,'(8(e16.6),i)') (x(ip)-1.)*r0, (y(ip)-1.)*r0,
     &           (z(ip)-1.)*r0, vx(ip)*v0/aexpn, vy(ip)*v0/aexpn,
     &           vz(ip)*v0/aexpn, pw(ip)*aM0*hubble, 0.0, ip
         endif
      enddo
c
      close(io)
      return
      END
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE Check_Contamination (xc, yc, zc, rmax, exitflag)
c--------------------------------------------------------------------
c     check the kind of DM particles in a halo
c
c     input in code units
      include 'a_tree.h'
      include 'a_control.h'
#     include 'a_def.h'
c
      integer i, ip
      real*8 xc, yc, zc, rmax, rmax2
      logical exitflag
      integer numparticles(nspecies)
      real*8 massparticles(nspecies)
      real*8 goodmass, badmass
      exitflag = .false.
      numparticles = 0
      massparticles = 0.
      goodmass = 0.
      badmass = 0.
c
      rmax2 = rmax**2
      do i = 1,nspecies
         do ip = nsp(i,1), nsp(i,2) !loop over all particles
            rpart = (x(ip)-xc)**2 + (y(ip)-yc)**2 + (z(ip)-zc)**2
            if (rpart .le. rmax2) then
               numparticles(i) = numparticles(i) + 1
               massparticles(i) = massparticles(i) + pw(ip)
            endif
         enddo
      enddo
c
      do i = 1,nspecies
         write(*,'(i,e14.6)') numparticles(i), massparticles(i)*aM0
      enddo

      goodmass = massparticles(1)
#ifdef STARFORM
      goodmass = goodmass + massparticles(nspecies)
      do i = 2,nspecies-1
         badmass = badmass + massparticles(i)
      enddo
#else
      do i = 2,nspecies
         badmass = badmass + massparticles(i)
      enddo
#endif
      write(*,'(2(a,e14.6,a))') ' bad = ',badmass*aM0,' Msun, ',
     &     ' good = ',goodmass*aM0,' Msun'
      if (badmass .gt. goodmass) then
         exitflag = .true.
      endif

      RETURN
      END
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE RegionTemperature(xc,yc,zc,rmax,Tgas)
      real*8 Tfact, T_code, T_code_sum
      Tfact = T_0 * (gamma - 1.) / (aexpn**2)
      T_code = 0.
      T_code_sum = 0.

      call Get_MaxLevelNow ()
      DO Level = MinLevel+1, MaxLevelNow  
         call Select_Cells ( Level , nLevel ) 
         do ic1 = 1 , nLevel
            icell = iSelect(ic1)
            do ic2 = 0 , 7
               ic = icell + ic2
               if ( iOctCh(ic) .eq. nil ) then !just leafs
                  T_code = hvar(8,ic)/hvar(1,ic)
                  T_code_sum = T_code_sum + T_code
               endif
            enddo
         enddo
      ENDDO

      Tgas = Tfact * T_code_sum ! K

      RETURN
      END
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE AllProfile(filename,xc,yc,zc,rmin,rmax,nrbin,
     &     Mtot,Mdm,Mgas,Mst,Tstars)
c-------------------------------------------------------------------- 
c     Make spherical profiles for gas properties and all components 
c     based on Output_Average_profile in ART_Analysis.F
c
c     input in code units.
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"

      character*256 filename 
      real*8 xc, yc, zc, rmin, rmax
      integer nrbin, nbins
      real*8 Mtot, Mdm, Mgas, Mst, Tstars
      real*8 gamma1, a3, a2, ai, a2i
      real*8 rho0C, rho0Cm, rlogmin, rlmax, drl,
     &     r0Ch, r0C, r0P, OmegaA, deltavir
      integer m, nDMparticles
      parameter ( npoints = 1000000 )
      double precision xx , yy, zz
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins)
      real dgbin(0:nmaxbins), dgibin(0:nmaxbins)
      real dstbin(0:nmaxbins) , dstibin(0:nmaxbins), nstbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), ndmp(0:nmaxbins)
#ifdef STARFORM
      real dpbin(nspecies-1,0:nmaxbins), dpibin(nspecies-1,0:nmaxbins),
     &     dpi(nspecies-1)
      integer npbin(nspecies-1,0:nmaxbins), npi(nspecies-1)
#endif
#ifndef STARFORM
      real dpbin(nspecies,0:nmaxbins), dpibin(nspecies,0:nmaxbins),
     &     dpi(nspecies)
      integer npbin(nspecies,0:nmaxbins), npi(nspecies)
#endif
      character*256 myformat
      integer nformat
      real dtibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins), entbin(0:nmaxbins)
      real gm(0:nmaxbins), gmc(0:nmaxbins)
      integer nbin(0:nmaxbins)

      integer io
      io = 42

      a3 = aexpn**3
      a2 = aexpn**2
      ai = 1.0 / aexpn
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      rho0Cm =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3 ! o n_H in cm^-3
      rho0C =  rho0 / a3 / 1.E+18     ! in M_sun/pc3 (proper units)
      r0Ch= r0*1000.0 ! in kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in kpc (Physical units) 
      OmegaA= Om0 / ( Om0 + Oml0*a3)
      deltavir = (18*pi**2+82*(OmegaA-1.) - 39*(OmegaA-1.)**2 )/OmegaA
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
#ifdef STARFORM
      nDMparticles=lsp(nspec-1)
#endif
#ifndef STARFORM
      nDMparticles=lsp(nspec)
#endif
!     initialize cell volumes
      do iLevel=MinLevel,MaxLevel
         CellVol(iLevel)= 2.0**(-3.0*iLevel) 
      enddo
! Initialize profiles: 
      do i = 0 , nmaxbins
         rbin(i) = 0. ! averaged radius inside shell
         dgbin(i) = 0. ! gas density
         ddmbin(i) = 0. ! DM density
         pbin(i) = 0. ! pressure
         tbin(i) = 0. ! temperature
         entbin(i) = 0. ! entropy
         nbin(i) = 0 ! number of virtual points inside shell 
         ndmp(i) = 0. ! number of DM particles inside shell
         dgibin(i)=0.  ! Gas mass inside radius r
         ddmibin(i) = 0. ! DM mass inside radius r
         do j = 1,nspecies
            npbin(j,i) = 0      ! n all dm
            dpbin(j,i) = 0.     ! dm density
            dpibin(j,i) = 0.    ! DM masses inside radius r
         enddo
         dtibin(i) = 0. ! Total mass inside radius r
         dstbin(i)=0. ! Stellar density 
         dstibin(i)=0. ! Stellar mass inside radius r
         gm(i)=0. ! Gas mass proved inside bin
         gmc(i)=0. ! Cold Gas mass proved inside bin
      enddo
!     Sample of random points inside a sphere of rmax:
      m = 28548535              !  seed
      call Get_MaxLevelNow ()
      DO i = 1 , npoints
         rr    = rmax * RANDd ( m )
         theta = pi2 * RANDd ( m ) - pi
         phi   = pi2 * RANDd ( m )
         xx = rr * cos(theta) * cos(phi) + xc
         yy = rr * cos(theta) * sin(phi) + yc
         zz = rr * sin(theta) + zc
!     Check periodic boundaries is not necesary.
         idcell = iFindCell ( MaxLevel , xx , yy , zz )
         pd= hvar(1,idcell)
         pp =  hvar(6,idcell)
         cw = pd * CellVol(iLv(idcell)) ! Cell mass

         if ( hvar(1,idcell) .gt. 0. ) then
            T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
            ptemp = T_0 * T_code * a2i
            pent  = S_0 * T_code / hvar(1,idcell)**gamma1
         else
            write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
         endif

         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            nbin(ibin)   = nbin(ibin) + 1
            rbin(ibin)   = rbin(ibin) + rp
            dgbin(ibin)  = dgbin(ibin) + pd
            pbin(ibin) = pbin(ibin) + pp*cw
            tbin(ibin) = tbin(ibin) + ptemp*cw
            entbin(ibin) = entbin(ibin) + pent*cw
            gm(ibin) = gm(ibin) + cw
            gmc(ibin) = gmc(ibin) + cw
         endif
      ENDDO
#ifdef STARFORM
      do j = 1,nspecies-1
#endif
#ifndef STARFORM
      do j = 1,nspecies
#endif
         do ip = nsp(j,1),nsp(j,2)
            xx = x(ip)
            yy = y(ip)
            zz = z(ip)
            rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
            ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
            if ( ibin .le. nmaxbins ) then
               npbin(j,ibin) = npbin(j,ibin) + 1 
               dpbin(j,ibin) = dpbin(j,ibin) + pw(ip)
            endif
         enddo
      enddo
      DO ip=1,nDMparticles ! DM
         xx = x(ip)
         yy = y (ip)
         zz = z (ip)
         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            ndmp(ibin) = ndmp(ibin) + 1
            ddmbin(ibin) = ddmbin(ibin) + pw(ip)
         endif
      ENDDO
#ifdef STARFORM
      DO ip=lsp(nspec-1)+1 , lsp(nspec)  ! Stars
         xx = x(ip)
         yy = y (ip)
         zz = z (ip)
         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            nstbin(ibin) = nstbin(ibin) + 1
            dstbin(ibin) = dstbin(ibin) + pw(ip)
         endif
      ENDDO
#endif
!     Normalize profiles: Mass weighted profiles for pressure, Temperature and entropy
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rnp = 1.d0 / gm(i)
            rbin(i) = rbin(i) / nbin(i) *r0P
            dgbin(i) = dgbin(i) / nbin(i)
            pbin(i) = pbin(i) * rnp 
            tbin(i) = tbin(i) *rnp
            entbin(i) = entbin(i) *rnp
         endif
      enddo
      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.                  ! DM mass inside radius
      dsi = 0.                   ! stellar mass inside radius
      dgi = 0.                   ! Gas mass inside radius
      dti = 0.                  ! Total mass inside radius
      dpi = 0.                  ! each particle mass inside radius
      DO i = 0 , nbins
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
            rl = 10.**(rlogmin + float(i-1)*drl)
         endif
         rr = 10.**(rlogmin + float(i)*drl)
         volr = 4.18879 * (rr**3 - rl**3)    
         if ( volr .gt. 0. ) then
            dmi = dmi + ddmbin(i) 
            dsi = dsi + dStbin(i)
            dgi = dgi + dgbin(i)*volr
            dti = dti + ddmbin(i) + dgbin(i)*volr + dStbin(i)

            ddmbin(i) = ddmbin(i)/ volr
            dStbin(i)  = dStbin(i) / volr
 
!            ddmibin(i) = dmi  / (4.18879 * rr**3) 
            ddmibin(i) = dmi 
            dStibin(i) = dsi  
            dgibin(i) = dgi  
            dtibin(i) = dti  / (4.18879 * rr**3)

#ifdef STARFORM
            do j = 1,nspecies-1
#endif
#ifndef STARFORM
            do j = 1,nspecies
#endif
               dpi(j) = dpi(j) + dpbin(j,i)
               dpbin(j,i) = dpbin(j,i) / volr
               dpibin(j,i) = dpi(j)
            enddo
         else
            write(*,*) 'bad volr =',volr,i, rl, rr
         endif
      ENDDO
c
c     calculate outer radius values
      Mtot = (dtibin(nbins)*4.18879*rho0/hubble**3)*
     &     (10.**(rlogmin + float(nbins)*drl)*r0)**3
      Mdm = ddmibin(nbins)*aM0
      Mgas = dgibin(nbins)*aM0
      Mst = dstibin(nbins)*aM0
      Tstars = 0.
#ifdef STARFORM      
      do is = 1, (nsp(nspecies,2) - nsp(nspecies,1) + 1)
         ip = is + nsp(nspecies,1)-1
         rp = sqrt((x(ip)-xc)**2+(y(ip)-yc)**2+(z(ip)-zc)**2)
         if (rp .le. rmax) then
            Tstars = Tstars + (pt(ip)-tbirth(is))
         endif
      enddo
      Tstars = Tstars/(nsp(nspecies,2)-nsp(nspecies,1)+1)
      Tstars = Tstars*t0
#endif
c-------------------------------------------------------------------- 
!     Write outputfile:
      open(io, file = TRIM(filename))
      write(io,10) (xc-1.)*r0  , (yc-1.)*r0 , (zc-1.)*r0    
      write(io,11) rmin*r0P, rmax*r0P , nrbin 
      write(io,12) rmax*r0P, Mtot
      write(io,14) aexpn , 1/aexpn -1.
 10   format ('# (xc,yc,zc)[Mpc h-1 comoving]= ', 3(F9.4,1x))
 11   format ('# rmin = ',g13.6,' rmax =',g13.6,
     &     ' nrbin =',i5)
 12   format ('# Rvir =',f8.1,
     &     ' [kpc], Mvir =',g13.6, ' [Msun]')
 14   format
     &('# Distances in kpc ,densities in cm-3 and M in M_Sun,
     &  for a= ',f6.3, ' and redshift= ', f6.2)
#ifdef STARFORM
      WRITE(myformat, '(3(a,i),a)') '(3(a,12x),a,10x,3(a,9x),a,10x,',
     &     (nspecies-1),'(14x),3(a,3x),',
     &     (nspecies-1),'(17x),a,4x,',
     &     (nspecies-1),'(4x),a)'
#endif
#ifndef STARFORM
      WRITE(myformat, '(3(a,i),a)') '(3(a,12x),a,10x,3(a,9x),a,10x,',
     &     (nspecies),'(14x),3(a,3x),',
     &     (nspecies),'(17x),a,4x,',
     &     (nspecies),'(4x),a)'
#endif
      write(io,myformat) '#      <r>',' r_l',' r_m',' r_r',
     &     'T gas (K)',' Rho_g ','Rho_St ','Rho_dm',
     &     'M_insideR_Gas','M_insideR_Stars',' M_insideR_DM ',
     &     'N_dm_Pibin',' M_insideR_tot'
      
#ifdef STARFORM
      nformat = 8+(nspecies-1)+3+(nspecies-1)
      WRITE(myformat, '(a,i,a,i,a)') '(',nformat,'(2x,g14.7),',
     &     (nspecies-1),'(2x,i7),1x,g14.7)'
#endif
#ifndef STARFORM
      nformat = 8+nspecies+3+nspecies
      WRITE(myformat, '(a,i,a,i,a)') '(',nformat,'(2x,g14.7),',
     &     nspecies,'(2x,i7),1x,g14.7)'
#endif

      NdmPibin=0.
      npi = 0                   ! number of particles inside radius
      do i = 0 , nbins
         NdmPibin = NdmPibin + ndmp(i)
#ifdef STARFORM
         do j = 1,nspecies-1
#endif
#ifndef STARFORM
         do j = 1,nspecies
#endif
            npi(j) = npi(j) + npbin(j,i)
         enddo

         if ( i .eq. 0 ) then 
            rl = 0.0  
            rmid = 0.5*10.**rlogmin 
         else
            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
         endif
         rr = 10.**(rlogmin + float(i)*drl) ! rightside
         if ( nbin(i) .gt. 0 ) then
            write(io,myformat)
     &           rbin(i), rl*r0P, rmid*r0P, rr*r0P, tbin(i), 
     &           dgbin(i)*rho0Cm, dStbin(i)*rho0Cm, ddmbin(i)*rho0Cm, 
#ifdef STARFORM
     &           (dpbin(j,i)*rho0Cm, j=1,nspecies-1),
#endif
#ifndef STARFORM
     &           (dpbin(j,i)*rho0Cm, j=1,nspecies),
#endif
     &           dgibin(i)*aM0, dstibin(i)*aM0, ddmibin(i)*aM0,
#ifdef STARFORM
     &           (dpibin(j,i)*aM0, j=1,nspecies-1),
     &           (npi(j), j=1,nspecies-1),
#endif
#ifndef STARFORM
     &           (dpibin(j,i)*aM0, j=1,nspecies),
     &           (npi(j), j=1,nspecies),
#endif
     &           dtibin(i)*4.18879*(rr*r0)**3*rho0/hubble**3
         endif
      enddo
c 17   format(13(2x,g14.7),1x,i7,1x,g14.7)
      close(io)
      return
      end
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE STARcenter (rmin, rmax,
     &     Xdmc , Ydmc , Zdmc, Vxdmc, Vydmc, Vzdmc)
c-------------------------------------------------------------------- 
c     Find center of the stellar distribution 
c     based on the iterative method in PMgalaxy.f
c
c  input in code units
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real*8 rmin, rmax, Xdmc , Ydmc , Zdmc, Vxdmc, Vydmc, Vzdmc
      real*8 x8, y8, z8,v8x, v8y, v8z
      Rselect=rmax
      Ncount = 0
 300  CONTINUE
      R2 = Rselect**2
      nn=0
      x8 =0.
      y8 = 0.
      z8 = 0.
      v8x=0.
      v8y =0.
      v8z =0.
      do ip= nsp(nspec,1),nsp(nspec,2)
         rp2 =(x(ip)-Xdmc)**2 + (y(ip)-Ydmc )**2 + (z(ip) - Zdmc)**2
         if (rp2.lt.R2) then
            nn = nn + 1
            x8 = x8 + x(ip)
            y8 = y8 + y(ip)
            z8 = z8 + z(ip)
            v8x = v8x + vx(ip)
            v8y = v8y + vy(ip)
            v8z = v8z + vz(ip)
         endif
      enddo
      If(nn.eq.0)write (*,*)' Error: no particles found. R=',Rselect
      Xdmc = x8/nn
      Ydmc = y8 / nn
      Zdmc = z8 / nn
      Vxdmc = v8x / nn
      Vydmc = v8y / nn
      Vzdmc = v8z / nn
c      write (*,12) Ncount,  Rselect*1000.0*aexpn*r0/hubble ,
c     &  (xdmc-1.)*1000.0*aexpn*r0/hubble,
c     &  (ydmc-1.)*1000.0*aexpn*r0/hubble,  
c     &   (zdmc-1.)*1000.0*aexpn*r0/hubble,
c     $ vxdmc* v0 , vydmc * v0 , vzdmc * v0,
c     & nn 
      Rselect = Rselect / 1.1
      Ncount = Ncount + 1
      write (*,*)  Rselect
      if ( Ncount.lt.1000.AND.Rselect.gt.rmin.AND.nn.gt.20) goto 300
 12   format (I3, 7(f10.3), I6)
      return
      end
c
c
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE angularM(xc,yc,zc,Vxc,Vyc,Vzc,isys,Rs,Tgas,aLx,aLy,aLz)
c-------------------------------------------------------------------- 
c
c     Calculate the total angular momentum inside a 
c     sphere Rs for a given center. 
c
c     isys specifies the system:
c     isys=1 --> DM particles
c     isys=2 --> stars
c     isys=3 --> gas below Tgas
c--------------------------------------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      real*8 xc, yc, zc, Vxc, Vyc, Vzc
      real*8 Rs, Tgas, aLx, aLy, aLz, aL
      real*8 Rs2, r0P, v0P
      real*8 sumx, sumy, sumz, aMassT, gamma1, a2, a2i, T_code
      integer isys, nn, nn2, nDMparticles
      integer ip, i1, icell, i2, idc
      real*8 xp, yp, zp, rp2, vpx, vpy, vpz
      real*8 Xpp, Ypp, Zpp, r2, pd, ptemp
      real   Posx, Posy, Posz
      real*8 Posx1, Posy1, Posz1
      real*8 vxpp, vypp, vzpp, aMassCell
c
      sumx=0.0
      sumy=0.0
      sumz=0.0
      aMassT=0.0
      nn=0
      nn2=0
      Rs2 = Rs**2
      r0P = r0 *1000.0 / hubble * aexpn ! in kpc (Physical units)
      v0P = v0 / aexpn          ! in km/s (proper units)
      gamma1 = gamma - 1.0
      a2 = aexpn**2
      a2i = 1.0 / a2
      IF (isys .eq. 1) then       ! DM particles:
         write (*,*) 'Angular momentum calculation for DM particles'
#ifdef STARFORM
      nDMparticles=lsp(nspecies-1)
#endif
#ifndef STARFORM
      nDMparticles=lsp(nspecies)
#endif
         DO ip=1,nDMparticles
            xp = x(ip) -xc
            yp = y(ip) -yc
            zp = z(ip) -zc
            rp2 = xp**2 + yp**2 + zp**2
            vpx =  vx(ip) - Vxc
            vpy =  vy(ip) - Vyc
            vpz =  vz(ip) - Vzc
            IF (rp2 .lt. Rs2) then
               nn= nn+1
               sumx = sumx + pw(ip)*( yp*vpz - zp*vpy)
               sumy = sumy + pw(ip)*( zp*vpx - xp*vpz) 
               sumz = sumz + pw(ip)*( xp*vpy - yp*vpx)
            ENDIF
         ENDDO
         write (*,'(a,f8.2,a,i)')  'The number of particles inside ',
     &               Rs*r0P,' proper kpc is ',nn
         aLx = sumx
         aLy = sumy
         aLz = sumz
      ELSEIF (isys .eq. 2) then   ! Stars
         write (*,*) 'Angular momentum calculation for star particles'
check         write (*,*)  nsp(nspec,1),nsp(nspec,2)
         DO ip=nsp(nspecies,1),nsp(nspecies,2)
            xp = x(ip) -xc
            yp = y(ip) - yc
            zp = z(ip) -zc
            rp2 = xp**2 + yp**2 + zp**2
            vpx =  vx(ip) - Vxc
            vpy =  vy(ip) - Vyc
            vpz =  vz(ip) - Vzc
            IF (rp2 .lt. Rs2) then
               nn = nn+1
               sumx = sumx + pw(ip)*( yp*vpz - zp*vpy)
               sumy = sumy + pw(ip)*( zp*vpx - xp*vpz) 
               sumz = sumz + pw(ip)*( xp*vpy - yp*vpx)
            ENDIF
         ENDDO
         write (*,'(a,f8.2,a,i)')  'The number of particles inside ',
     &               Rs*r0P, ' proper kpc is ',nn
         aLx = sumx
         aLy = sumy
         aLz = sumz
      ELSEIF (isys .eq. 3) then   ! Gas below some Tgas
         write (*,'(a,a,g12.4,a)') 'Angular momentum calculation',
     &     ' for gas with T <',Tgas,' K'
c
         call Get_MaxLevelNow ()
         DO Level =   MinLevel+1, MaxLevelNow 
            CellVoLevel= 2.0**(-3.0*Level)
            nLevel = iNOLL(Level)
            call Select_Cells ( Level , nLevel ) 
            do i1 = 1 , nLevel
               icell = iSelect(i1)
               do i2 = 0 , 7
                  idc = icell + i2
                  if ( iOctCh(idc) .eq. nil ) then ! only leaves
                     call Ps ( idc, Posx, Posy, Posz )
                     Posx1 = Posx
                     Posy1 = Posy
                     Posz1 = Posz
                     call PeriodicBoundaryD2(xc,yc,zc,Posx1,Posy1,Posz1)
                     Xpp = Posx1-xc
                     Ypp = Posy1-yc
                     Zpp = Posz1-zc
c                     call PeriodicBoundary(xc,yc,zc,Posx,Posy,Posz)
c                     Xpp = Posx-xc
c                     Ypp = Posy-yc
c                     Zpp = Posz-zc
                     r2 = (Xpp)**2+(Ypp)**2+(Zpp)**2
c                     if( Level .eq. MaxLevelNow) then
c                     if (i2 .eq. 0) then
c                        write(*,'(6(f8.4,1x))') Posx1,xc,Posy1,yc,Posz1,zc
c                        write(*,'(3(f10.5,1x))') Xpp, Ypp, Zpp
c                        write(*,*) sqrt(Posx1**2+Posy1**2+Posz1**2)*r0P,
c     &                   sqrt(r2)*r0P
c                        write(*,*) '---'
c                     endif
c                     endif
                     if (r2 .lt. Rs2) then ! cell center inside rs
                        nn2=nn2+1
                        pd= hvar(1,idc)
                        if (pd .gt. 0.) then
                           T_code= gamma1 * hvar(8,idc)/hvar(1,idc)
                           ptemp = T_0 * T_code * a2i
                        else
                           write(*,*)'bad gas density in cell',idc,pd
                        endif
                        IF (ptemp .lt. Tgas) then
                           nn=nn+1
c
                           vxpp= hvar(3,idc)/pd - Vxc
                           vypp= hvar(4,idc)/pd - Vyc
                           vzpp= hvar(5,idc)/pd - Vzc
c
                           aMassCell=pd*CellVoLevel
                           aLx =(Ypp*Vzpp-Zpp*Vypp)*aMassCell
                           aLy =(Zpp*Vxpp-Xpp*Vzpp)*aMassCell
                           aLz =(Xpp*Vypp-Ypp*Vxpp)*aMassCell
c
                           sumx = sumx + aLx
                           sumy = sumy + aLy
                           sumz = sumz + aLz
c
                           aMassT=aMassT+aMassCell
                        ENDIF
                     endif
                  endif
               enddo
            enddo
         ENDDO         
         write (*,'(a,f8.2,a,i)')  'The number of cells inside ',Rs*r0P,
     &                ' proper kpc is ',nn2
         write (*,'(a,f8.2,a,i)')  'The number of cells also with T < ',
     &                Tgas,' is ',nn
         aLx = sumx
         aLy = sumy
         aLz = sumz                           
c
         write (*,'(a,f12.5)') 'The specific angular momentum is',
     &            sqrt(aLx**2+aLy**2+aLz**2)/aMassT *r0P*v0P
      ENDIF
      END
C-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE NewAxes(aLx,aLy,aLz,AxX,AxY,AxZ)
c-------------------------------------------------------------------- 
c
c     Find new system of coordinates 
c     in which aM L=(aLx,aLy,aLz) is parallel to V3
c
c     based on PMgalaxy.f
c
c     input in code units
c-------------------------------------------------------------------- 
      real*8 aLx, aLy, aLz, AxX(3) , AxY(3), AxZ(3)
      real*8 L, Lx, Ly, Lz, Ld
      real*8 ex1, ey1, ez1, ex2, ey2, ez2, ex3, ey3, ez3
      real*8 prod1, prod2, prod3, a1, a2, a3
c
c
! Find new system of coordinates: ex,ey,ez
      L = sqrt( aLx**2 +aLy**2 + aLz**2)
      Lx = aLx/L
      Ly = aLy/L
      Lz = aLz/L
      Ld = sqrt(Lx**2 +Lz**2)
      ex1 =  Lz/Ld
      ey1 = 0.
      ez1 = -Lx/Ld
      ex2 = -Lx*abs(Ly)/Ld
      ey2 =  Ld
      ez2 = -Lz*abs(Ly)/Ld
c         If(Ly.lt.0.) ey2 = -Ld
c         If(Ly.lt.0.) then
c            ey2 = -Ld
c            ex1=-ex1
c            ey1=-ey1
c            ez1=-ez1
c         endif
      IF (Ly.lt.0.) then
         ex2=-ex2
         ez2=-ez2
      ENDIF
c            
      ex3 = Lx
      ey3 = Ly
      ez3 = Lz
c                     Check if the vectors are (1) unit (2) orthogonal
      prod1 = ex1*Lx +ey1*Ly +ez1*Lz
      prod2 = ex2*Lx +ey2*Ly +ez2*Lz
      prod3 = ex1*ex2 +ey1*ey2 +ez1*ez2
      a1    = Lx**2 +Ly**2 +Lz**2
      a2    = ex1**2 +ey1**2 +ez1**2
      a3    = ex2**2 +ey2**2 +ez2**2
c      write (*,*)  'New system of coordinates: Vector 1, 2 ,3' 
c      write (*,*) ' Vector 1=',ex1,ey1,ez1 
c      write (*,*) ' Vector 2=',ex2,ey2,ez2 
c      write (*,*) ' Vector 3=',ex3,ey3,ez3 
c      write (*,*) ' Vector lengths=', a1,a2,a3,' must be 1'
c      write (*,*) ' Vector products=', prod1,prod2,prod3,' must be 0'
      if( (abs(a1-1.).gt.1.e-2).or.(abs(a2-1.).gt.1.e-2).or.
     &    (abs(a2-1.).gt.1.e-2).or.(abs(prod1).gt.1.5e-1).or.
     &    (abs(prod2).gt.1.5e-1).or.(abs(prod3).gt.1.5e-1) ) then
         write (*,*) 'ERROR: broken at NewAxes conversion'
         stop
      endif
      AxX(1) = ex1
      AxX(2) = ey1
      AxX(3) = ez1
      AxY(1) = ex2
      AxY(2) = ey2
      AxY(3) = ez2
      AxZ(1) = ex3
      AxZ(2) = ey3
      AxZ(3) = ez3
      return
      end
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE Write_Stars_ASCII(filename,xc,yc,zc,Vxc,Vyc,Vzc,rmax)
     &     
c    spherical distribution with current mass and metallicity
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"

      character*256 filename
      real*8 xc, yc, zc, Vxc, Vyc, Vzc, rmax
      real*8 Xph, Yph, Zph, Vxph, Vyph, Vzph
      real*8 aMassTotal, aMassTotal0 
      real*8 ageS, ageU, age
      real*8 xfact, vfact, tfact, rmax2
      integer is, istar, Nstars
      integer io
      io = 43
c
      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      tfact = t0*aexpn**2
      rmax2 = rmax**2
c
      open(io, file=TRIM(filename))
      write(io,300) t*tfact, aexpn
 300  format ('# time=', E20.6, ' yrs , a= ', F10.6)
      write(io,301) (xc-1.)*r0 ,(yc-1.)*r0 ,(zc-1.)*r0
 301  format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
      write(io,'(a,a)') '#X, Y, Z (kpc), Vx, Vy, Vz (km/s), ',
     &     'm0, m (Msun), age (Gyr), SNII, SIa (metal mass fraction)'
c
      Nstars = 0
      aMassTotal = 0.
      aMassTotal0 = 0.
      do is = nsp(nspecies,1), nsp(nspecies,2) ! stars
         istar = is - nsp(nspecies,1) + 1
         Xph = x(is) - xc
         Yph = y(is) - yc
         Zph = z(is) - zc
         Vxph = vx(is) - Vxc
         Vyph = vy(is) - Vyc
         Vzph = vz(is) - Vzc
c
         IF ((Xph**2 + Yph**2 + Zph**2) .lt. rmax2) THEN 
            ageS = ageU - age(dble( tbirth(istar))) 
            Nstars = Nstars+1 
            aMassTotal = aMassTotal +  pw(is) 
            aMassTotal0 = aMassTotal0 +  pw0(istar) 
            write(io,334) istar, Xph*xfact, Yph*xfact, Zph*xfact,
     &           Vxph*vfact, Vyph*vfact, Vzph*vfact,
     &           pw0(istar)*aM0, pw(is)*aM0, ageS, 
     &           zstII(istar), zstIa(istar)
 334        format(I7, 2x, 11(e14.6,1x))
         ENDIF
      enddo
      close(io)
c
      write(*,'(a,f12.4,a,i)') ' N stars inside ',rmax*xfact,
     &     ' kpc :',Nstars                   
      write(*,'(a,e12.4,a)') ' Total mass = ', aMassTotal*aM0, 
     &     ' Msun'
      write(*,'(a,e12.4,a)') ' Total initial mass = ', aMassTotal0*aM0,
     &     ' Msun'
      end
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE Write_Stars_ASCII2(fprofname,xc,yc,zc,Vxc,Vyc,Vzc,
     &     rmax,tAge,AxX,AxY,AxZ,iopt)
c    iopt 0 = spherical distribution with current mass
c    iopt 1 = spherical distribution with initial mass
c    iopt 2 = cylindrical distribution with current mass
c    iopt 3 - no option
c    iopt 4 = spherical distribution with current mass and metallicity
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real*8 xc,yc,zc,Vxc,Vyc,Vzc,rmax,tAge
      real AxX(3), AxY(3), AxZ(3)            ! New axes for cylindrical profiles
      real*8 Xph,Yph,Zph, Vxph, Vyph, Vzph
      character*256 fprofname 
      real*8 aMassTotal0 
      real*8 aMassTotal, aMcum, aMcumlow, aMcumhigh !AM distribution
      parameter ( nmaxbins = 1000 )
      parameter (aLzlow=0.1e3)
      real aMtotal(nmaxbins),aMLzlow(nmaxbins), aMLzhigh(nmaxbins)
      real*8 ageS, ageU, age,tAgeGyr
      integer io
      io = 43

      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      aL0 = aM0*vfact*xfact 
      aL0s= vfact*xfact
      tfact = t0*aexpn**2
      rmax2= rmax**2
! Age:
      ageU = age(t)             ! in Gyr
      tAgeGyr =  tAge / 1.e9    ! in Gyr  

!AM distribution:
      drl=rmax/100.
      aMassTotal=0.0
      do i = 1 , nmaxbins
         aMtotal(i)=0.0
         aMLzlow(i)=0.0
         aMLzhigh(i)=0.0
      enddo
c
      open(io, file= TRIM(fprofname))
      write(io,300) t*tfact, aexpn
      write(io,301) (xc-1.)*r0 ,(yc-1.)*r0 ,(zc-1.)*r0
      if (iopt.eq.4) then
         write(io,303)
      else
         write(io,302)
      endif
      Nstars=0
      aMassTotal0=0.0
      do is=nsp(nspecies,1),nsp(nspecies,2) ! stars
         istar =is - nsp(nspecies,1) + 1
         Xph =  x(is)-xc
         Yph =  y(is)-yc
         Zph =  z(is)-zc

         Vxph = vx(is) - Vxc
         Vyph = vy(is) - Vyc
         Vzph = vz(is) - Vzc

         IF (IOPT.EQ.4) THEN    ! Spherical distribution with current mass and metallicites 
            IF ( Xph**2 + Yph**2 + Zph**2 .lt. rmax2) THEN   
!               if ( (t-tbirth(istar))*tfact.lt.tAge) then  
               ageS = ageU - age(dble( tbirth(istar))) 
               if (ageS.le.tAgeGyr) then
                  Nstars=Nstars+1  
                  aMassTotal=aMassTotal +  pw(is)  
                  aMassTotal0=aMassTotal0 +  pw0(istar) 
                  write(io,334) istar,
     &                 (Xph)*xfact , 
     &                 (Yph)*xfact , 
     &                 (Zph)*xfact , 
     &                 Vxph*vfact, Vyph*vfact, Vzph*vfact, 
     &                 pw(is)*aM0 ,                                                          
!     &                 pw0(istar)*aM0 ,   
     &                 ageS,
     &                 zstII(istar),zstIa(istar)
                   endif   
            ENDIF

         ELSEIF (IOPT.EQ.0) THEN ! Spherical distribution with current mass 
            IF ( Xph**2 + Yph**2 + Zph**2 .lt. rmax2) THEN  
!               if ( (t-tbirth(istar))*tfact.lt.tAge) then 
               ageS = ageU - age(dble( tbirth(istar)))  
               if (ageS.le.tAgeGyr) then 
                  Nstars=Nstars+1  
                  aMassTotal=aMassTotal +  pw(is) 
                  aMassTotal0=aMassTotal0 +  pw0(istar)  
                  write(io,334) istar, 
     &                 (Xph)*xfact , 
     &                 (Yph)*xfact ,  
     &                 (Zph)*xfact , 
c    &                 vx(is)*vfact,vy(is)*vfact,vz(is)*vfact, 
     &                 Vxph*vfact, Vyph*vfact, Vzph*vfact, 
     &                 pw(is)*aM0 , 
c     &                 pw0(istar)*aM0 ,                                               
c     &                 (t-tbirth(istar))*tfact 
     &                 ageS 
                   endif 
            ENDIF 

         ELSEIF (IOPT.EQ.1) THEN    ! Spherical distribution with initial mass
            IF ( Xph**2 + Yph**2 + Zph**2 .lt. rmax2) THEN
!               if ( (t-tbirth(istar))*tfact.lt.tAge) then
               ageS = ageU - age(dble( tbirth(istar)))
               if (ageS.le.tAgeGyr) then
                  Nstars=Nstars+1
                  aMassTotal=aMassTotal +  pw(is)
                  aMassTotal0=aMassTotal0 +  pw0(istar)      
                  write(io,334) istar,
     &                 (Xph)*xfact , 
     &                 (Yph)*xfact ,
     &                 (Zph)*xfact ,
c    &                 vx(is)*vfact,vy(is)*vfact,vz(is)*vfact,
     &                 Vxph*vfact, Vyph*vfact, Vzph*vfact,
c     &                 pw(is)*aM0 , 
     &                 pw0(istar)*aM0 , 
c     &                 (t-tbirth(istar))*tfact
     &                 ageS
                   endif
            ENDIF
         ELSEIF (IOPT.EQ.2) THEN ! Cilindrical distribution
            rp3D2 = xph**2 + yph**2 + zph**2 
            Xpp =xph*AxX(1) +yph*AxX(2) +zph*AxX(3)
            Ypp =xph*AxY(1) +yph*AxY(2) +zph*AxY(3)
            Zpp =xph*AxZ(1) +yph*AxZ(2) +zph*AxZ(3)
!
            Vxpp=Vxph*AxX(1) +Vyph*AxX(2) +Vzph*AxX(3)
            Vypp=Vxph*AxY(1) +Vyph*AxY(2) +Vzph*AxY(3)
            Vzpp=Vxph*AxZ(1) +Vyph*AxZ(2) +Vzph*AxZ(3)
!
            aLx =(Ypp*Vzpp-Zpp*Vypp)!*pw(is)
            aLy =(Zpp*Vxpp-Xpp*Vzpp)!*pw(is)
            aLz =(Xpp*Vypp-Ypp*Vxpp)!*pw(is)

!
            rp3D2p = xpp**2 + ypp**2 + zpp**2 
            if (abs( rp3D2 -  rp3D2p).gt.0.01) then
           write (*,*) 'error in the rotation of coordinates. Rp3D2= ',
     &              rp3D2, rp3D2p  
           write (*,*)  'it will not stop'
!     STOP
            endif
            rp2D=  xpp**2 + ypp**2
            IF (rp3D2p.lt.rmax2) THEN ! AM DISTRIBUTION
!               if ( (t-tbirth(istar))*tfact.lt.tAge) then
               ageS = ageU - age(dble( tbirth(istar)))
               if (ageS.le.tAgeGyr) then
                  Nstars=Nstars+1
                  write(io,334) istar,
     &                 Xpp*xfact , 
     &                 Ypp*xfact ,
     &                 Zpp*xfact ,
!     &                 vx(is)*vfact,vy(is)*vfact,vz(is)*vfact,
     &                 vxpp*vfact, vypp*vfact, vzpp*vfact,
     &                 pw(is)*aM0 , 
     &                 aLx*aL0s, aLy*aL0s, aLz*aL0s,
!     &                 (t-tbirth(istar))*tfact
     &                 ageS
! AM distribution:
                  ibin = int(sqrt(rp3D2p)/drl)+1
                  aMtotal(ibin)=aMtotal(ibin)+pw(is)
                  if (aLz*aL0s.lt.aLzlow) then
                     aMLzlow(ibin)=aMLzlow(ibin)+pw(is)
                  else
                     aMLzhigh(ibin)=aMLzhigh(ibin)+pw(is)
                  endif
                  aMassTotal = aMassTotal +pw(is)
                  aMassTotal0=aMassTotal0 +  pw0(istar)
               endif
            ENDIF
         ENDIF
      enddo
      write(io,'(a,f8.2,a,i)') ' N stars inside ',rmax*xfact,' kpc :',
     &           Nstars
      write(io,*) 'Total mass:',  aMassTotal*aM0
      write(io,*) 'Total initial mass:', aMassTotal0*aM0
      write(io,340) tAge
c
      write(*,*)' N stars inside ',rmax*xfact,' kpc :',Nstars                   
      write(*,*) ' Total mass:',  aMassTotal*aM0                               
      write(*,*) ' Total initial mass:', aMassTotal0*aM0                       
      write(*,340) tAge
      write(*,*) ' SFR= ', aMassTotal0*aM0/tAge
c
 300  format ( 'time=', E20.6, ' yrs , a= ', F10.6)
 301  format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
 302  format('X Y Z [kpc],Vx  Vy Vz [ km/s] M0[Msun] age[Gyr]')
 303  format('X Y Z [kpc],Vx  Vy Vz [ km/s] M0[Msun] age[Gyr] 
     &   Mass fraction in SNII,SIa metals')
 334  format(I7 , 2x, 11(e14.6,1x))
 340  format ( ' younger than ', E8.1, ' yrs')
c
! AM distribution
      aMcum=0.0
      aMcumlow=0.0
      aMcumhigh=0.0
      do i=1,nmaxbins
         rr = float(i)*drl
         aMcum=aMcum + aMtotal(i)
         aMcumlow= aMcumlow + aMLzlow(i)
         aMcumhigh=aMcumhigh+ aMLzhigh(i)
         write (io,341) i, rr*xfact, 
     &      aMcumlow/aMassTotal, aMcumhigh/aMassTotal,aMcum/aMassTotal 
      enddo
 341  format(I4, 4(g14.3))
      close(io) 
      end
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE ProfilesSp(filename,xc,yc,zc,Tgas,rmin,rmax,nrbin,
     &                      rvirHalo,aMvirHalo)
c-------------------------------------------------------------------- 
c     Make spherical profiles for gas properties and all components 
c     based on Output_Average_profile in ART_Analysis.F
c
c     input in code units.
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"

      real*8 xc, yc, zc, Tgas, rmin, rmax, rvirHalo, aMvirHalo
      integer nrbin
      character*256 filename 
      real*8 gamma1, a3, a2, ai, a2i
      real*8 rho0C, rho0Cm, rlogmin, rlmax, drl,
     &     r0Ch, r0C, r0P, OmegaA, deltavir
      integer m, nDMparticles
      parameter ( npoints = 1000000 )
      double precision xx , yy, zz
      parameter ( nmaxbins = 1000 )
      real rbin(0:nmaxbins)
      real dgbin(0:nmaxbins), dgibin(0:nmaxbins)
      real dstbin(0:nmaxbins) , dstibin(0:nmaxbins), nstbin(0:nmaxbins)
      real ddmibin(0:nmaxbins), ddmbin(0:nmaxbins), ndmp(0:nmaxbins)
      real dtibin(0:nmaxbins),  bfibin(0:nmaxbins)
      real tbin(0:nmaxbins),pbin(0:nmaxbins), entbin(0:nmaxbins)
      real gm(0:nmaxbins), gmc(0:nmaxbins)
      integer nbin(0:nmaxbins)

      integer io
      io = 44

      a3 = aexpn**3
      a2 = aexpn**2
      ai = 1.0 / aexpn
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      rho0Cm =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3 ! o n_H in cm^-3
      rho0C =  rho0 / a3 / 1.E+18     ! in M_sun/pc3 (proper units)
      r0Ch= r0*1000.0 ! in kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in kpc (Physical units) 
      OmegaA= Om0 / ( Om0 + Oml0*a3)
!      deltavir = 119           ! with respect to the matter mean density (code units) at z>>1
!      deltavir = 340. ! same but at z=0
      deltavir = (18*pi**2+82*(OmegaA-1.) - 39*(OmegaA-1.)**2 )/OmegaA
c
c      write (*,*)  'Overdensity 500 critical'
c      deltavir=deltavir*5.
c   
      rlogmin = log10(rmin)
      rlmax = log10(rmax)
      drl = (rlmax - rlogmin) / float(nrbin) 
      nbins = int((rlmax-rlogmin)/drl) + 1
#ifdef STARFORM
      nDMparticles=lsp(nspec-1)
#endif
#ifndef STARFORM
      nDMparticles=lsp(nspec)
#endif
!     initialize cell volumes
      do iLevel=MinLevel,MaxLevel
         CellVol(iLevel)= 2.0**(-3.0*iLevel) 
      enddo
! Initialize profiles: 
      do i = 0 , nmaxbins
         rbin(i) = 0. ! averaged radius inside shell
         dgbin(i) = 0. ! gas density
         ddmbin(i) = 0. ! DM density
         pbin(i) = 0. ! pressure
         tbin(i) = 0. ! temperature
         entbin(i) = 0. ! entropy
         nbin(i) = 0 ! number of virtual points inside shell 
         ndmp(i) = 0. ! number of DM particles inside shell
         bfibin(i)=0            ! baryonic fraction
         dgibin(i)=0.  ! Gas mass inside radius r
         ddmibin(i) = 0. ! DM mass inside radius r
         dtibin(i) = 0. ! Total mass inside radius r
         dstbin(i)=0. ! Stellar density 
         dstibin(i)=0. ! Stellar mass inside radius r
         gm(i)=0. ! Gas mass proved inside bin
         gmc(i)=0. ! Cold Gas mass proved inside bin
      enddo
!     Sample of random points inside a sphere of rmax:
      m = 28548535              !  seed
      call Get_MaxLevelNow ()
      DO i = 1 , npoints
         rr    = rmax * RANDd ( m )
         theta = pi2 * RANDd ( m ) - pi
         phi   = pi2 * RANDd ( m )
         xx = rr * cos(theta) * cos(phi) + xc
         yy = rr * cos(theta) * sin(phi) + yc
         zz = rr * sin(theta) + zc
!     Check periodic boundaries is not necesary.
         idcell = iFindCell ( MaxLevel , xx , yy , zz )
         pd= hvar(1,idcell)
         pp =  hvar(6,idcell)
         cw = pd * CellVol(iLv(idcell)) ! Cell mass

         if ( hvar(1,idcell) .gt. 0. ) then
            T_code = gamma1 * hvar(8,idcell)/hvar(1,idcell)
            ptemp = T_0 * T_code * a2i
            pent  = S_0 * T_code / hvar(1,idcell)**gamma1
         else
            write(*,*) 'bad gas density in cell',idcell,hvar(1,idcell)
         endif

         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            nbin(ibin)   = nbin(ibin) + 1
            rbin(ibin)   = rbin(ibin) + rp
            dgbin(ibin)  = dgbin(ibin) + pd
            pbin(ibin) = pbin(ibin) + pp*cw
            tbin(ibin) = tbin(ibin) + ptemp*cw
            entbin(ibin) = entbin(ibin) + pent*cw
            gm(ibin) = gm(ibin) + cw
            IF (ptemp.lt.Tgas) then
               gmc(ibin) = gmc(ibin) + cw
            ENDIF
         endif
      ENDDO
      DO ip=1,nDMparticles ! DM
         xx = x(ip)
         yy = y (ip)
         zz = z (ip)
         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            ndmp(ibin) = ndmp(ibin) + 1
            ddmbin(ibin) = ddmbin(ibin) + pw(ip)
         endif
      ENDDO
#ifdef STARFORM
      DO ip=lsp(nspec-1)+1 , lsp(nspec)  ! Stars
         xx = x(ip)
         yy = y (ip)
         zz = z (ip)
         rp = sqrt((xx-xc)**2+(yy-yc)**2+(zz-zc)**2)
         ibin = max ( int((log10(rp) - rlogmin)/drl) + 1 , 0 )
         if ( ibin .le. nmaxbins ) then
            nstbin(ibin) = nstbin(ibin) + 1
            dstbin(ibin) = dstbin(ibin) + pw(ip)
         endif
      ENDDO
#endif
!     Normalize profiles: Mass weighted profiles for pressure, Temperature and entropy
      do i = 0 , nbins
         if ( nbin(i) .gt. 0 ) then
            rnp = 1.d0 / gm(i)
            rbin(i) = rbin(i) / nbin(i) *r0P
            dgbin(i) = dgbin(i) / nbin(i)
            pbin(i) = pbin(i) * rnp 
            tbin(i) = tbin(i) *rnp
            entbin(i) = entbin(i) *rnp
         endif
      enddo
      rhodmavei = 1.0 / (1.0 - Omb0 / Om0) ! average DM density
      dmi = 0.                  ! DM mass inside radius
      dsi =0.                   ! stellar mass inside radius
      dgi =0.                   ! Gas mass inside radius
      dti = 0.                  ! Total mass inside radius
      DO i = 0 , nbins
         if ( i .eq. 0 ) then 
            rl = 0.0
         else
            rl = 10.**(rlogmin + float(i-1)*drl)
         endif
         rr = 10.**(rlogmin + float(i)*drl)
         volr = 4.18879 * (rr**3 - rl**3)    
         if ( volr .gt. 0. ) then
            dmi = dmi + ddmbin(i) 
            dsi = dsi + dStbin(i)
            dgi = dgi + dgbin(i)*volr
            dti = dti + ddmbin(i) + dgbin(i)*volr + dStbin(i)

            ddmbin(i) = ddmbin(i)/ volr
            dStbin(i)  = dStbin(i) / volr
 
!            ddmibin(i) = dmi  / (4.18879 * rr**3) 
            ddmibin(i) = dmi 
            dStibin(i) = dsi  
            dgibin(i) = dgi  
            dtibin(i) = dti  / (4.18879 * rr**3)

           bfibin(i) = (dgi*(gmc(i)/gm(i))+dsi) / dti  
!	bfibin(i) = dsi / dti   ! stellar to total mass ratio       
            if ( i .gt. 0 .and.
     &           dtibin(i) .lt. deltavir .and.
     &           dtibin(i-1) .ge. deltavir ) then 
               ivir = i
               rvir = (deltavir * (rr - rl) + rl*dtibin(i) - 
     &              rr*dtibin(i-1)) / (dtibin(i) - dtibin(i-1))
c     rvir = rr 
            endif
            if ( i .gt. 0 .and.
     &           ddmibin(i) .lt. deltavir .and.
     &           ddmibin(i-1) .ge. deltavir ) then
            ivirDM = i
            rvirDM = (deltavir * (rr - rl) + rl*ddmibin(i) -
     &           rr*ddmibin(i-1)) / (ddmibin(i) - ddmibin(i-1))
c     rvir = rr 
         endif  
      else
         write(*,*) 'bad volr =',volr,i, rl, rr
      endif
      ENDDO
      aMvir = deltavir * 
     &     4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      aMvirDM = deltavir *
     &     4.18879 * (rvirDM*r0)**3 * rho0 / hubble**2
c-------------------------------------------------------------------- 
      rvirHalo=rvir*r0P
      aMvirHalo=aMvir/hubble 
c--------------------------------------------------------------------
!     Write outputfile:
      open(io, file = TRIM(filename))
c      write(io,12) ivirDM, rvirDM*r0P, aMvirDM/hubble, deltavir
      write(io,10) (xc-1.)*r0  , (yc-1.)*r0 , (zc-1.)*r0    
      write(io,11) rmin*r0P, rmax*r0P , nrbin 
      write(io,12) ivir, rvir*r0P, aMvir/hubble, deltavir
      write(io,13) c
      write(io,14) aexpn , 1/aexpn -1. , Tgas
      write(io,15)
 10   format ('# (xc,yc,zc)[Mpc h-1 comoving]= ', 3(F9.4,1x))
 11   format ('# rmin = ',g13.6,' rmax =',g13.6,
     &     ' nrbin =',i5)
 12   format ('# ivir = ',i4,' Rvir =',f8.1,
     &     ' [kpc], Mvir =',g13.6,
     &     ' [Msun] for Dvir(mean)=',f7.2)
 13   format ('# NFW fit: concentration= ', f3.0)
 14   format
     &('# Distances in kpc ,densities in cm-3 and M in M_Sun,
     &  for a= ',f6.3, ' and redshift= ', f6.2, 
     &' Baryonic fraction of cold gas', E8.1, ' K ')
 15   format ('# <r>',14x,' r_l',14x,' r_m',14x,' r_r',14x,
     &     ' Rho_g ',3x, 'Rho_St ',8x, 'Rho_dm',10x,
     &     'M_insideR_Gas',3x,  'M_insideR_Stars',3x,
     &     ' M_insideR_DM ',3x,' M_insideR_tot',4x,
     &     'T_g',10x,'S_g',8x,
     &     'N_dm_Pibin',2x,'Baryonic_F  [Uni]' )
!     16   format ('#  h^-1 Mpc  overdensities w.r.t. mean dens.',
!     &        '             K  keV cm^2 ')
      NdmPibin=0.
      do i = 0 , nbins
         NdmPibin =NdmPibin + ndmp(i)
         if ( i .eq. 0 ) then 
            rl = 0.0  
            rmid = 0.5*10.**rlogmin 
         else
            rl = 10.**(rlogmin + float(i-1)*drl) ! leftside
            rmid = 10.**(rlogmin + float(i)*drl - 0.5*drl) ! middle
         endif
         rr = 10.**(rlogmin + float(i)*drl) ! rightside
         if ( nbin(i) .gt. 0 ) then
            write(io,17)
     &           rbin(i), rl*r0P  , 
     &         rmid*r0P , rr*r0P  , 
     &           dgbin(i)*rho0Cm ,dStbin(i)*rho0Cm, ddmbin(i)*rho0Cm, 
     &           dgibin(i)*aM0   ,
     &           dstibin(i)*aM0    ,
     &           ddmibin(i)*aM0 ,
     &           dtibin(i)*4.18879*(rr*r0)**3*rho0/hubble**3 ,  
     &           tbin(i),entbin(i),int(NdmPibin), bfibin(i)/(Omb0/Om0)
!     &           tbin(i),entbin(i),int(nstbin(i)), bfibin(i)/(Omb0/Om0) 
!    dgibin(i)*4.18879*(rr*r0)**3*rho0/hubble**3
         endif
      enddo
 17   format(13(2x,g14.7),1x,i7,1x,g14.7)
      close(io)
      write(*,*) 'done writing spherical profiles.'
      return
      end
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE Write_DM_binary_zoom(filename,xc,yc,zc,
     &     Vxc,Vyc,Vzc,rmax)
      !writes out a binary file with all the DM particles
      ! the file is centered on the selected halo
      ! data output: id, x, y, z, vx, vy, vz, and mass 
      !  for each particle in physical units: kpc km/s and Msun
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
c
      character*256 filename 
      real*8 xc,yc,zc,Vxc,Vyc,Vzc,rmax
      real*8 Xph,Yph,Zph,Vxph,Vyph,Vzph
      real*8 xfact, vfact, tfact, rmax2, aMassTotal
      integer is, Ndm
      integer io
      io = 45
c
      write(*,*) ''
      write(*,'(a,a)') ' writing ', TRIM(filename)
c
      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      tfact = t0*aexpn**2
c
      open(io, file=TRIM(filename), form='unformatted')
      Ndm=0
      aMassTotal=0.
      do is=nsp(1,1), nsp(nspecies-1,2) ! all DM species
         Xph = x(is)-xc
         Yph = y(is)-yc
         Zph = z(is)-zc
         Vxph = vx(is)-Vxc
         Vyph = vy(is)-Vyc
         Vzph = vz(is)-Vzc
         IF((ABS(Xph).lt.rmax).AND.(ABS(Yph).lt.rmax).AND.
     &        (ABS(Zph).lt.rmax))THEN
            Ndm=Ndm+1
            aMassTotal = aMassTotal + pw(is)
            write(io) is, Xph*xfact, Yph*xfact, Zph*xfact,
     &           Vxph*vfact, Vyph*vfact, Vzph*vfact,
     &           pw(is)*aM0
         ENDIF
      enddo
      close(io) 
c
      write(*,'(a,a,f12.4,a,i)') ' N DM particles inside a box with ',
     &     'sides of ',2.*rmax*xfact,' kpc :',Ndm
      write(*,'(a,e12.4,a)') ' Tot mass inside = ',aMassTotal*aM0,
     &     ' Msun'
      end
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE Write_Stars_binary_zoom(filename,xc,yc,zc,
     &     Vxc,Vyc,Vzc,rmax)
      !writes out a binary file with all the stellar particles
      ! the file is centered on the selected halo
      ! data output: id, x, y, z, vx, vy, vz, initial mass, current mass, age, SNII, SNIa
      !  for each particle in physical units: kpc, km/s, Msun and metallicity mass fraction
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"

      character*256 filename 
      real*8 xc,yc,zc,Vxc,Vyc,Vzc,rmax
      real*8 Xph,Yph,Zph, Vxph, Vyph, Vzph
      real*8 aMassTotal0, aMassTotal
      real*8 ageS, ageU, age
      integer io
      io = 46
c
      write(*,*) ''
      write(*,'(a,a)') ' writing ', TRIM(filename)
c
      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      tfact = t0*aexpn**2
      ageU = age(t)             ! in Gyr
c
      open(io, file=TRIM(filename), form = 'unformatted')
      Nstars = 0
      aMassTotal = 0.
      aMassTotal0 = 0.
      do is = nsp(nspecies,1),nsp(nspecies,2) ! stars
         istar = is - nsp(nspecies,1) + 1
         Xph = x(is)-xc
         Yph = y(is)-yc
         Zph = z(is)-zc
         Vxph = vx(is)-Vxc
         Vyph = vy(is)-Vyc
         Vzph = vz(is)-Vzc
         IF ((ABS(Xph).lt.rmax).AND.(ABS(Yph).lt.rmax).AND.
     &    (ABS(Zph).lt.rmax)) THEN
            ageS = ageU - age(dble( tbirth(istar)))
            Nstars = Nstars+1
            aMassTotal = aMassTotal + pw(is)
            aMassTotal0 = aMassTotal0 + pw0(istar)      
            write(io) istar, Xph*xfact, Yph*xfact, Zph*xfact,
     &           Vxph*vfact, Vyph*vfact, Vzph*vfact,
     &           pw0(istar)*aM0, pw(is)*aM0, sngl(ageS),
     &           zstII(istar), zstIa(istar)
         ENDIF
      enddo
      close(io)
c
      write(*,'(a,a,f12.4,a,i)') ' N star particles inside a box ',
     &     'with sides of ',2.*rmax*xfact,' kpc :',Nstars
      write(*,'(a,e12.4,a)') ' Total mass = ',  aMassTotal*aM0, ' Msun'  
      write(*,'(a,e12.4,a)') ' Total initial mass = ', aMassTotal0*aM0,  
     &     ' Msun'
      end
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE Write_Gas_binary_zoom(filename,xc,yc,zc,
     &     Vxc,Vyc,Vzc,rmax)
      ! Write out the gas information as a binary file
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
c
      character*256 filename
      real*8 xc,yc,zc,Vxc,Vyc,Vzc,rmax
      real Posx, Posy, Posz
      real*8 Posxc, Posyc, Poszc, aMassTotal
      integer Ncells, i1, icell, i2, idc
      real*8 T_code, ptemp, xfact, gamma1, a2, a2i, ai,
     &     r0Pk, v0P, rho0Cm, Zsolar
      integer io
      io = 47
c
      write(*,*) ''
      write(*,'(a,a)') ' writing ', TRIM(filename)
c
      xfact = 1000.0*aexpn*r0/hubble
      gamma1 = gamma - 1.0 
      a2 = aexpn**2
      a2i = 1.0 / a2
      ai = 1.0 / aexpn
      r0Pk= r0*1000.0 / hubble * aexpn ! in kpc (Physical units)
      v0P = v0 / aexpn          ! in km/s (proper units)
      rho0Cm =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3 ! o n_H in cm^-3
      Zsolar = 0.0199 
c
      call Get_MaxLevelNow ()
      Ncells = 0
      aMassTotal = 0.
      open(io, file=TRIM(filename), form='unformatted')
      DO Level = MaxLevelNow, MinLevel+1,-1 
         nLevel = iNOLL(Level)
         call Select_Cells ( Level , nLevel ) 
         do i1 = 1 , nLevel
            icell = iSelect(i1)
            do i2 = 0 , 7
               idc = icell + i2
               if ( iOctCh(idc) .eq. nil ) then ! only leaves
                  call  Ps(idc, Posx, Posy, Posz)
                  Posxc = Posx - xc
                  Posyc = Posy - yc
                  Poszc = Posz - zc
c
                  IF ((abs(Posxc).le.rmax).AND.(abs(Posyc).le.rmax).AND.
     &                 (abs(Poszc).le.rmax)) THEN
                     rho = hvar(1,idc)
                     vx0 = hvar(3,idc) / rho - Vxc
                     vy0 = hvar(4,idc) / rho - Vyc
                     vz0 = hvar(5,idc) / rho - Vzc
                     T_code = gamma1 * hvar(8,idc)/rho
                     ptemp = T_0 * T_code * a2i
                     ZII = hvar(izII,idc)/rho ! /Zsolar
                     ZIa = hvar(izIa,idc)/rho ! /Zsolar
c     
                     write(io) 2.**(-Level)*r0Pk*1000.,
     &                    Posxc*r0Pk ,Posyc*r0Pk, Poszc*r0Pk,  
     &                    vx0*v0P,vy0*v0P,vz0*v0P, rho*rho0Cm, ptemp,
     &                    ZII, ZIa
                     Ncells = Ncells+1
                     aMassTotal = aMassTotal + rho*CellSize(Level)**3
                  ENDIF
               endif
            enddo
         enddo
      ENDDO
      Level=0
      DO i1 = 1 , ncell0
         iCh1 = iOctCh(i1)
         if ( iOctCh(i1) .eq. nil ) then ! only leaves
c        
            call  Ps ( i1 , Posx,Posy,Posz )
            Posxc = Posx - xc
            Posyc = Posy - yc
            Poszc = Posz - zc
            IF ((abs(Posxc).le.rmax).AND.(abs(Posyc).le.rmax).AND.
     &           (abs(Poszc).le.rmax)) then
               rho = hvar(1,i1)
               vx0 = hvar(3,i1) / rho - Vxc
               vy0 = hvar(4,i1) / rho - Vyc
               vz0 = hvar(5,i1) / rho - Vzc
               T_code = gamma1 * hvar(8,i1)/rho
               ptemp = T_0 * T_code * a2i
               ZII = hvar(izII,i1)/rho ! /Zsolar
               ZIa = hvar(izIa,i1)/rho ! /Zsolar
c
               write(io)  2.**(-Level)*r0Pk*1000.,
     &              Posxc*r0Pk ,Posyc*r0Pk, Poszc*r0Pk,  
     &              vx0*v0P, vy0*v0P, vz0*v0P, rho*rho0Cm, ptemp,
     &              ZII, ZIa
               Ncells = Ncells+1
               aMassTotal = aMassTotal + rho*CellSize(Level)**3
            ENDIF
         endif
      ENDDO
      close(io)
c
      write(*,'(a,f12.4,a,i)') ' N cells inside a box with sides of ',
     &     2.*rmax*xfact,' kpc :',Ncells
      write(*,'(a,e12.4,a)') ' Total mass = ', aMassTotal*aM0,' Msun'
      END
c--------------------------------------------------------------------
c
c--------------------------------------------------------------------
      SUBROUTINE Write_Gas_ASCII_zoom(filename,xc,yc,zc,
     &     Vxc,Vyc,Vzc,rmax)
      ! Write out the gas information as an ASCII file
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
c
      character*256 filename
      real*8 xc,yc,zc,Vxc,Vyc,Vzc,rmax
      real Posx, Posy, Posz
      real*8 Posxc, Posyc, Poszc, aMassTotal
      integer Ncells, i1, icell, i2, idc
      real*8 T_code, ptemp, xfact, gamma1, a2, a2i, ai,
     &     r0Pk, v0P, rho0Cm, Zsolar
      integer io
      io = 48
c
      write(*,*) ''
      write(*,'(a,a)') ' writing ', TRIM(filename)
c
      xfact = 1000.0*aexpn*r0/hubble
      gamma1 = gamma - 1.0
      a2 = aexpn**2
      a2i = 1.0 / a2
      ai = 1.0 / aexpn
      r0Pk= r0*1000.0 / hubble * aexpn ! in kpc (Physical units)
      v0P = v0 / aexpn          ! in km/s (proper units)
      rho0Cm =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3 ! o n_H in cm^-3
      Zsolar = 0.0199
c
      call Get_MaxLevelNow ()
      Ncells = 0
      aMassTotal = 0.
      open(io, file=TRIM(filename))
      write(io, '(5x,9(a,15x),a2,14x,a2)')
     &     '1','2','3','4','5','6','7','8','9','10','11'
      write(io, '(a10,11x,3(a1,15x),2(a2,14x),
     &      a2,11x,a7,8x,a11,3x,a14,3x,a14)')
     &     'cell size','x','y','z','vx','vy','vz','density',
     &     'temperature','SNII mass frac','SNIa mass frac'

      DO Level = MaxLevelNow, MinLevel+1,-1
         nLevel = iNOLL(Level)
         call Select_Cells ( Level , nLevel )
         do i1 = 1 , nLevel
            icell = iSelect(i1)
            do i2 = 0 , 7
               idc = icell + i2
               if ( iOctCh(idc) .eq. nil ) then ! only leaves
                  call  Ps(idc, Posx, Posy, Posz)
                  Posxc = Posx - xc
                  Posyc = Posy - yc
                  Poszc = Posz - zc
c
                  IF ((abs(Posxc).le.rmax).AND.(abs(Posyc).le.rmax).AND.
     &                 (abs(Poszc).le.rmax)) THEN
                     rho = hvar(1,idc)
                     vx0 = hvar(3,idc) / rho - Vxc
                     vy0 = hvar(4,idc) / rho - Vyc
                     vz0 = hvar(5,idc) / rho - Vzc
                     T_code = gamma1 * hvar(8,idc)/rho
                     ptemp = T_0 * T_code * a2i
                     ZII = hvar(izII,idc)/rho ! /Zsolar
                     ZIa = hvar(izIa,idc)/rho ! /Zsolar
c
                     write(io, '(es11.4,10es16.4)')
     &                    2.**(-Level)*r0Pk*1000.,
     &                    Posxc*r0Pk, Posyc*r0Pk, Poszc*r0Pk,
     &                    vx0*v0P, vy0*v0P,vz0*v0P, rho*rho0Cm, ptemp,
     &                    ZII, ZIa
                     Ncells = Ncells+1
                     aMassTotal = aMassTotal + rho*CellSize(Level)**3
                  ENDIF
               endif
            enddo
         enddo
      ENDDO
      Level=0
      DO i1 = 1 , ncell0
         iCh1 = iOctCh(i1)
         if ( iOctCh(i1) .eq. nil ) then ! only leaves
c
            call  Ps ( i1 , Posx,Posy,Posz )
            Posxc = Posx - xc
            Posyc = Posy - yc
            Poszc = Posz - zc
            IF ((abs(Posxc).le.rmax).AND.(abs(Posyc).le.rmax).AND.
     &           (abs(Poszc).le.rmax)) then
               rho = hvar(1,i1)
                      vx0 = hvar(3,i1) / rho - Vxc
               vy0 = hvar(4,i1) / rho - Vyc
               vz0 = hvar(5,i1) / rho - Vzc
               T_code = gamma1 * hvar(8,i1)/rho
               ptemp = T_0 * T_code * a2i
               ZII = hvar(izII,i1)/rho ! /Zsolar
               ZIa = hvar(izIa,i1)/rho ! /Zsolar
c
               write(io, '(es11.4,10es16.4)')
     &              2.**(-Level)*r0Pk*1000.,
     &              Posxc*r0Pk ,Posyc*r0Pk, Poszc*r0Pk,
     &              vx0*v0P, vy0*v0P, vz0*v0P, rho*rho0Cm, ptemp,
     &              ZII, ZIa
               Ncells = Ncells+1
               aMassTotal = aMassTotal + rho*CellSize(Level)**3
            ENDIF
         endif
      ENDDO
      close(io)
c
      write(*,'(a,f12.4,a,i)') ' N cells inside a box with sides of ',
     &     2.*rmax*xfact,' kpc :',Ncells
      write(*,'(a,e12.4,a)') ' Total mass = ', aMassTotal*aM0,' Msun'
      END
c-------------------------------------------------------------------- 
c
c-------------------------------------------------------------------- 
      SUBROUTINE Extract_BOX3(filename,xc,yc,zc,Vxc,Vyc,Vzc,QSOBox)
c-------------------------------------------------------------------- 
c
c     Extract info from all cells inside a QSOBox centered on galaxy
c
c     iflagHI=3 --> x,y,z, vx,vy,vz, density , temperature
c
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      character*256 filename
      real*8 rs, rsd, rhogl, T_g, Z_met, xHII, xHI, 
     &     coolrate, cool, heat, gamma1
      integer io
      io = 49

!      parameter (QSOBox=28) !300 ) ! 400 ) ! kpc
      QSOBoxcu = QSOBox/ ( r0 * 1000.0 / hubble * aexpn) ! in code units
      QSOBoxcu2= QSOBoxcu/2.0
      r0Pk= r0*1000.0 / hubble * aexpn ! in kpc (Physical units)
! For HI densities:
      d0 = dble(zero)
      d1 = dble(one)
      rs = max ( 1.d0 / aexpn - 1.0 , 0.d0 ) 
      gamma1 = gamma - 1.0 
      Zsolar = 0.0199 
      a2 = aexpn**2
      a2i = 1.0 / a2
      ai = 1.0 / aexpn
      fact_nH = log10(1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3) ! factor converting gas density in code units to n_H in cm^-3
!For physical units
      rho0Cm =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3 ! o n_H in cm^-3
      v0P = v0 / aexpn          ! in km/s (proper units)
c-------------------------------------------------------------------- 
      call Get_MaxLevelNow ()

c
      Ncells=0
c
      open(io, file = TRIM(filename), form = 'unformatted')

      DO Level = MaxLevelNow, MinLevel+1,-1 
         nLevel = iNOLL(Level)
         call Select_Cells ( Level , nLevel ) 
         do i1 = 1 , nLevel
            icell = iSelect(i1)
            do i2 = 0 , 7
               idc = icell + i2
               if ( iOctCh(idc) .eq. nil ) then ! only leaves 
                  
                  call  Ps ( idc , Posx,Posy,Posz )
                  Posxc = Posx - xc
                  Posyc = Posy - yc
                  Poszc = Posz - zc
            IF (abs(Posxc).le.QSOBoxcu2.AND.abs(Posyc).le.QSOBoxcu2.AND.
     &                 abs(Poszc).le.QSOBoxcu2) then

                   rho = hvar(1,idc)
                   vx0 =    hvar(3,idc) / rho - Vxc
                   vy0 =    hvar(4,idc) / rho - Vyc
                   vz0 =    hvar(5,idc) / rho - Vzc
                   T_code = gamma1 * hvar(8,idc)/rho
                   ptemp = T_0 * T_code * a2i
c     
               write(io) 2.**(-Level)*r0Pk*1000.,
     &              Posxc*r0Pk ,Posyc*r0Pk, Poszc*r0Pk,  
     &              vx0*v0P,vy0*v0P,vz0*v0P, rho*rho0Cm, ptemp
               Ncells=Ncells+1
            ENDIF
         endif
      enddo
      enddo
      ENDDO
      Level=0
      DO i1 = 1 , ncell0
         iCh1 = iOctCh(i1)
         if ( iOctCh(i1) .eq. nil ) then ! only leaves
c
            call  Ps ( i1 , Posx,Posy,Posz )
            Posxc = Posx - xc
            Posyc = Posy - yc
            Poszc = Posz - zc
            IF (abs(Posxc).le.QSOBoxcu2.AND.abs(Posyc).le.QSOBoxcu2.AND.
     &           abs(Poszc).le.QSOBoxcu2) then

                  rho = hvar(1,i1)
                  vx0 =    hvar(3,i1) / rho - Vxc
                  vy0 =    hvar(4,i1) / rho - Vyc
                  vz0 =    hvar(5,i1) / rho - Vzc
                  T_code = gamma1 * hvar(8,i1)/rho
                  ptemp = T_0 * T_code * a2i
c     
                  write(io)   2.**(-Level)*r0Pk*1000.,
     &                 Posxc*r0Pk ,Posyc*r0Pk, Poszc*r0Pk,  
     &                 vx0*v0P,vy0*v0P,vz0*v0P, rho*rho0Cm, ptemp
                  Ncells=Ncells+1

            ENDIF
         endif
      ENDDO
      write (*,*)  'Number of cells= ', Ncells
      close(io)
      END
c-------------------------------------------------------------------- 
      SUBROUTINE Extract_BOX6(filename,xc,yc,zc,Vxc,Vyc,Vzc,QSOBox)
c-------------------------------------------------------------------- 
c
c     Extract info from all cells inside a QSOBox centered on galaxy
c
c     iflagHI=6 --> + metallicites
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      character*256 filename
      real*8 rs, rsd, rhogl, T_g, Z_met, xHII, xHI, 
     & coolrate, cool, heat, gamma1
      integer io
      io = 49

!      parameter (QSOBox=28) !300 ) ! 400 ) ! kpc

      write(*,*) 'Extract_BOX6'

      QSOBoxcu = QSOBox/ ( r0 * 1000.0 / hubble * aexpn) ! in code units
      QSOBoxcu2= QSOBoxcu/2.0
      r0Pk= r0*1000.0 / hubble * aexpn ! in kpc (Physical units)
! For HI densities:
      d0 = dble(zero)
      d1 = dble(one)
      rs = max ( 1.d0 / aexpn - 1.0 , 0.d0 ) 
      gamma1 = gamma - 1.0 
      Zsolar = 0.0199 
      a2 = aexpn**2
      a2i = 1.0 / a2
      ai = 1.0 / aexpn
      fact_nH = log10(1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3) ! factor converting gas density in code units to n_H in cm^-3
!For physical units
      rho0Cm =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3 ! o n_H in cm^-3
      v0P = v0 / aexpn          ! in km/s (proper units)
c-------------------------------------------------------------------- 
      call Get_MaxLevelNow ()

c
      Ncells=0  
c 
      open(io, file = TRIM(filename), form = 'unformatted')

      DO Level = MaxLevelNow, MinLevel+1,-1 
         nLevel = iNOLL(Level)
         call Select_Cells ( Level , nLevel ) 
         do i1 = 1 , nLevel
            icell = iSelect(i1)
            do i2 = 0 , 7
               idc = icell + i2
               if ( iOctCh(idc) .eq. nil ) then ! only leaves
                  
                  call  Ps ( idc , Posx,Posy,Posz )
                  Posxc = Posx - xc
                  Posyc = Posy - yc
                  Poszc = Posz - zc
            IF (abs(Posxc).le.QSOBoxcu2.AND.abs(Posyc).le.QSOBoxcu2.AND.
     &                 abs(Poszc).le.QSOBoxcu2) then

                  rho = hvar(1,idc)
                  vx0 =    hvar(3,idc) / rho - Vxc
                  vy0 =    hvar(4,idc) / rho - Vyc
                  vz0 =    hvar(5,idc) / rho - Vzc
                  T_code = gamma1 * hvar(8,idc)/rho
                  ptemp = T_0 * T_code * a2i
c     
                  ZII =hvar(izII,idc)/rho ! /Zsolar
                  ZIa =hvar(izIa,idc)/rho ! /Zsolar
c 
                  write(io) 2.**(-Level)*r0Pk*1000.,
     &                 Posxc*r0Pk ,Posyc*r0Pk, Poszc*r0Pk,  
     &                 vx0*v0P,vy0*v0P,vz0*v0P, rho*rho0Cm, ptemp,
     &                 ZII, ZIa
                  Ncells=Ncells+1
               ENDIF
            endif
         enddo
      enddo
      ENDDO
      Level=0
      DO i1 = 1 , ncell0
         iCh1 = iOctCh(i1)
         if ( iOctCh(i1) .eq. nil ) then ! only leaves
c
            call  Ps ( i1 , Posx,Posy,Posz )
            Posxc = Posx - xc
            Posyc = Posy - yc
            Poszc = Posz - zc
            IF (abs(Posxc).le.QSOBoxcu2.AND.abs(Posyc).le.QSOBoxcu2.AND.
     &           abs(Poszc).le.QSOBoxcu2) then

                  rho = hvar(1,i1)
                  vx0 =    hvar(3,i1) / rho - Vxc
                  vy0 =    hvar(4,i1) / rho - Vyc
                  vz0 =    hvar(5,i1) / rho - Vzc
                  T_code = gamma1 * hvar(8,i1)/rho
                  ptemp = T_0 * T_code * a2i
c     
                  ZII =hvar(izII,i1)/rho ! /Zsolar
                  ZIa =hvar(izIa,i1)/rho ! /Zsolar
c
                  write(io)   2.**(-Level)*r0Pk*1000.,
     &                 Posxc*r0Pk ,Posyc*r0Pk, Poszc*r0Pk,  
     &                 vx0*v0P,vy0*v0P,vz0*v0P, rho*rho0Cm, ptemp,
     &                 ZII, ZIa
                  Ncells=Ncells+1
               ENDIF
            endif
         ENDDO
      write (*,*)  'Number of cells= ', Ncells
      close(io)
      END
c
c     Previously from analysis_ART.F
c
C---------------------
c	PROFILES
C---------------------
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
      real*8 xc, yc, zc, rmax
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
c     -----------------------------------------------------------
      subroutine Output_Average_Profiles ( xc , yc , zc , rvir, 
     &                            rmin , rmax , nrbin , ibintype ,  
     &                            deltavir, ifindcenter , 
     &                            lpointsout, l_force_center,
     &                            filename )
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
c     character*256 filename - full path and file name to output profile to 
c                               it must end with a space
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real xc , yc , zc , rmin , rmax , deltavir, gamma1
      integer nrbin, ibintype, ifindcenter 
      logical lpointsout, l_force_center 
      character*256 filename 
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

      integer io
      io = 49

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
         open(io, file = 'points.dat')
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
            write(io,*) rp*r0, pd(ii), ptemp(ii),pent(ii)
          endif
        endif
      enddo
      if ( lpointsout ) then       
         close(io)
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
      open(io, file = TRIM(filename))
c
      aMvir = deltavir * 
     &        4.18879 * (rvir*r0)**3 * rho0 / hubble**2 
      write(io,10) xc , yc , zc  
      write(io,11) rmin , rmax , nrbin , ibintype 
      write(io,12) ivir, rvir*r0, aMvir, deltavir
      write(io,13) ifindcenter, l_force_center
      write(io,14)
      write(io,15)
      write(io,16)
 10   format ('# (xc,yc,zc) = ', 3(g14.7,1x))
 11   format ('# rmin = ',g13.6,' rmax =',g13.6,
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
          write(io,17)
     &         rbin(i), rl*r0, rmid, rr*r0, 
     &         dgbin(i),ddmbin(i),ddmibin(i), dtibin(i), 
     &         tbin(i),entbin(i),int(ndmp(i))
	endif
      enddo
 17   format(10(2x,g14.7),1x,i7)
      close(io)
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
c
c     -----------------------------------------------------------
      subroutine PeriodicBoundaryD2 ( xc, yc, zc, Posx, Posy, Posz)
c     Take care of periodic boundary conditions for analysis
c     -----------------------------------------------------------
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 xc, yc, zc, Posx, Posy, Posz
      real*8 diff_x, diff_y, diff_z, corr_x, corr_y, corr_z
c     
      diff_x = xc - Posx
      diff_y = yc - Posy
      diff_z = zc - Posz
      corr_x = 0.0
      corr_y = 0.0
      corr_z = 0.0
      
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

