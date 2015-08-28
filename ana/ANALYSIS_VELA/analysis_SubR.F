c   Subroutines for analysis of cosmological simulations ( recursive version)
c
c     by Daniel Ceverino, Jerusalem (2010)
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

      real Xdmc , Ydmc , Zdmc, Vxdmc, Vydmc, Vzdmc
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
c-------------------------------------------------------------------- 
      SUBROUTINE angularM(vcx,vcy,vcz,xc,yc,zc,isys,Rs,Tgas,aLx,aLy,aLz)
c-------------------------------------------------------------------- 
c
c     Calculate the total angular momentum inside a 
c     sphere Rs for a given center. 
c     isys specifies the system:
c
c     isys=1 --> DM particles
c     isys=2 --> stars
c     isys=3 --> gas bellow T=10^4 K
c--------------------------------------------------------------------
      include 'a_tree.h'
      include 'a_control.h'
#     include "a_def.h"
      double precision sumx, sumy, sumz,aMassT
      sumx=0.0
      sumy=0.0
      sumz=0.0
      aMassT=0.0
      nn=0
      Rs2 = Rs**2
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
      v0P = v0 / aexpn          ! in km/s (proper units)
      gamma1 = gamma - 1.0
      a2 = aexpn**2
      a2i = 1.0 / a2
      IF (isys.eq.1) then       ! DM particles:
#ifdef STARFORM
      nDMparticles=lsp(nspecies-1)
#endif
#ifndef STARFORM
      nDMparticles=lsp(nspecies)
#endif
         DO ip=1,nDMparticles
            xp = x(ip) -xc
            yp = y(ip) - yc
            zp = z(ip) -zc
            rp2 = xp**2 + yp**2 + zp**2
            vpx =  vx(ip) - vcx
            vpy =  vy(ip) - vcy
            vpz =  vz(ip) - vcz
            IF ( rp2.lt.Rs2) then
               nn= nn+1
               sumx = sumx + pw(ip)*( yp*vpz - zp*vpy)
               sumy = sumy + pw(ip)*( zp*vpx - xp*vpz) 
               sumz = sumz + pw(ip)*( xp*vpy - yp*vpx)
            ENDIF
         ENDDO
         write (*,*)  'The number of particles inside ',Rs*r0P,' is ',nn
         aLx = sumx
         aLy = sumy
         aLz = sumz
      ELSEIF (isys.eq.2) then   ! Stars
check         write (*,*)  nsp(nspec,1),nsp(nspec,2)
         DO ip=nsp(nspecies,1),nsp(nspecies,2)
            xp = x(ip) -xc
            yp = y(ip) - yc
            zp = z(ip) -zc
            rp2 = xp**2 + yp**2 + zp**2
            vpx =  vx(ip) - vcx
            vpy =  vy(ip) - vcy
            vpz =  vz(ip) - vcz
            IF ( rp2.lt.Rs2) then
               nn = nn+1
               sumx = sumx + pw(ip)*( yp*vpz - zp*vpy)
               sumy = sumy + pw(ip)*( zp*vpx - xp*vpz) 
               sumz = sumz + pw(ip)*( xp*vpy - yp*vpx)
            ENDIF
         ENDDO
         write (*,*)  'The number of particles inside ',Rs*r0P,' is ',nn
            aLx = sumx
            aLy = sumy
            aLz = sumz
      ELSEIF (isys.eq.3) then   ! Gas bellow 10^4 K.
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
                     call  Ps ( idc , Posx,Posy,Posz )
                     Xpp=Posx-Xc
                     Ypp=Posy-Yc
                     Zpp=Posz-Zc
                     r2=(Xpp)**2+(Ypp)**2+(Zpp)**2
                     if (r2.lt.rs2) then ! cell inside rs
                        pd= hvar(1,idc)
                        if ( pd .gt. 0. ) then
                           T_code= gamma1 * hvar(8,idc)/hvar(1,idc)
                           ptemp = T_0 * T_code * a2i
                        else
                           write(*,*)'bad gas density in cell',idc,pd
                        endif
                        IF (ptemp.lt.Tgas) then
                           nn=nn+1
c
                           vxpp= hvar(3,idc)/pd - vcx
                           vypp= hvar(4,idc)/pd - vcy
                           vzpp= hvar(5,idc)/pd - vcz
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
         write (*,*)  'The number of cells inside ',Rs*r0P,' is ',nn
         aLx = sumx
         aLy = sumy
         aLz = sumz                           
c
         write (*,*)  'The specific angular momentum of cold gas is'
         write (*,*)  sqrt(aLx**2+aLy**2+aLz**2)/aMassT *r0P*v0P
      ENDIF
      END
C-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE NewAxes(aLx,aLy,aLz,v1,v2,v3)
c-------------------------------------------------------------------- 
c
c     Find new system of coordinates 
c     in which aM L=(aLx,aLy,aLz) is parallel to V3
c
c     based on PMgalaxy.f
c
c     input in code units
c-------------------------------------------------------------------- 
      real Lx, Ly, Lz, L , Ld
      real v1(3) , v2(3), v3(3)

! Find new system of coordinates: ex,ey,ez
      L = sqrt( aLx**2 +aLy**2 + aLz**2)
      Lx = aLx/L
      Ly = aLy/L
      Lz = aLz/L
      Ld =sqrt(Lx**2 +Lz**2)
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
            
      ex3 =Lx
      ey3 =Ly
      ez3 =Lz
c                     Check if the vectors are (1) unit (2) orthogonal
      prod1 = ex1*Lx +ey1*Ly +ez1*Lz
      prod2 = ex2*Lx +ey2*Ly +ez2*Lz
      prod3 = ex1*ex2 +ey1*ey2 +ez1*ez2
      a1       = Lx**2 +Ly**2 +Lz**2
      a2       = ex1**2 +ey1**2 +ez1**2
      a3       = ex2**2 +ey2**2 +ez2**2
      write (*,*)  'New system of coordinates: Vector 1, 2 ,3' 
      write (*,*) ' Vector 1=',ex1,ey1,ez1 
      write (*,*) ' Vector 2=',ex2,ey2,ez2 
      write (*,*) ' Vector 3=',ex3,ey3,ez3 
      write (*,*) ' Vector lengths=', a1,a2,a3,' must be 1'
      write (*,*) ' Vector products=', prod1,prod2,prod3,' must be 0'
      if(abs(a1-1.).gt.1.e-2)stop
      if(abs(a2-1.).gt.1.e-2)stop
      if(abs(a2-1.).gt.1.e-2)stop
      if(abs(prod1).gt.1.5e-1)stop
      if(abs(prod2).gt.1.5e-1)stop
      if(abs(prod3).gt.1.5e-1)stop
      v1(1) = ex1
      v1(2) = ey1
      v1(3) = ez1
      v2(1) = ex2
      v2(2) = ey2
      v2(3) = ez2
      v3(1) = ex3
      v3(2) = ey3
      v3(3) = ez3
      return
      end
c-------------------------------------------------------------------- 
      SUBROUTINE Write_Stars_ASCII(fprofname, Vxc, Vyc, Vzc,
     &     xc,yc,zc, rmax,tage,ZDmax,AxX,AxY,AxZ,iopt)
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmax,tage,ZDmax,Vxc, Vyc, Vzc
      real*8 Xph,Yph,Zph, Vxph, Vyph, Vzph
      character*256 fprofname 
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles
      real*8 aMassTotal0
!AM distribution:
      real*8 aMassTotal, aMcum, aMcumlow, aMcumhigh
      parameter ( nmaxbins = 1000 )
      parameter (aLzlow=0.1e3)
      real aMtotal(nmaxbins),aMLzlow(nmaxbins), aMLzhigh(nmaxbins)
      real*8 ageS, ageU, age,tageGyr

      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      aL0 = aM0*vfact*xfact 
      aL0s= vfact*xfact
      tfact = t0*aexpn**2
      rmax2= rmax**2
! Age:
      ageU = age(t)             ! in Gyr
      tageGyr =  tage / 1.e9    ! in Gyr  

!AM distribution:
      drl=rmax/100.
      aMassTotal=0.0
      do i = 1 , nmaxbins
         aMtotal(i)=0.0
         aMLzlow(i)=0.0
         aMLzhigh(i)=0.0
      enddo

      nfn  = index(fprofname, ' ') - 1
      open ( 33, file=fprofname(1:nfn) )

      write(33,300) t*tfact, aexpn
      write(33,301) (xc-1.)*r0 ,(yc-1.)*r0 ,(zc-1.)*r0
      if (iopt.eq.4) then
         write(33,303)
      else
         write(33,302)
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
!               if ( (t-tbirth(istar))*tfact.lt.tage) then                                       
               ageS = ageU - age(dble( tbirth(istar)))                                           
               if (ageS.le.tageGyr) then                                                         
                  Nstars=Nstars+1                                                                
                  aMassTotal=aMassTotal +  pw(is)                                                
                  aMassTotal0=aMassTotal0 +  pw0(istar)                                          
                  write(33,334) istar,                                                           
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
!               if ( (t-tbirth(istar))*tfact.lt.tage) then                                       
               ageS = ageU - age(dble( tbirth(istar)))                                           
               if (ageS.le.tageGyr) then                                                         
                  Nstars=Nstars+1                                                                
                  aMassTotal=aMassTotal +  pw(is)                                                
                  aMassTotal0=aMassTotal0 +  pw0(istar)                                          
                  write(33,334) istar,                                                           
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
!               if ( (t-tbirth(istar))*tfact.lt.tage) then
               ageS = ageU - age(dble( tbirth(istar)))
               if (ageS.le.tageGyr) then
                  Nstars=Nstars+1
                  aMassTotal=aMassTotal +  pw(is)
                  aMassTotal0=aMassTotal0 +  pw0(istar)      
                  write(33,334) istar,
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
!            IF (rp2D.lt.rmax2.AND.ABS(ZPP).LT.ZDmax) then
            IF (rp3D2p.lt.rmax2) THEN ! AM DISTRIBUTION
!               if ( (t-tbirth(istar))*tfact.lt.tage) then
               ageS = ageU - age(dble( tbirth(istar)))
               if (ageS.le.tageGyr) then
                  Nstars=Nstars+1
                  write(33,334) istar,
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
      write(33,*)"Nstars inside ",rmax*xfact," Kpc :",Nstars
      write(33,*) "Total mass:",  aMassTotal*aM0
      write(33,*) "Total initial mass:", aMassTotal0*aM0
      write(33,340) tage
c
      write(*,*)"Nstars inside ",rmax*xfact," Kpc :",Nstars                   
      write(*,*) "Total mass:",  aMassTotal*aM0                               
      write(*,*) "Total initial mass:", aMassTotal0*aM0                       
      write(*,340) tage
      write(*,*) 'SFR= ', aMassTotal0*aM0/tage
c
 300  format ( 'time=', E20.6, ' yrs , a= ', F10.6)
 301  format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
 302  format("X Y Z [Kpc],Vx  Vy Vz [ km/s] M0[Msun] age[Gyr]")
 303  format("X Y Z [Kpc],Vx  Vy Vz [ km/s] M0[Msun] age[Gyr] 
     &   Mass fraction in SNII,SIa metals")
 334  format(I7 , 2x, 11(e14.6,1x))
 340  format ( ' younger than ', E6.1, ' yrs')
      close ( 33 ) 
! AM distribution
      aMcum=0.0
      aMcumlow=0.0
      aMcumhigh=0.0
      do i=1,nmaxbins
         rr = float(i)*drl
         aMcum=aMcum + aMtotal(i)
         aMcumlow= aMcumlow + aMLzlow(i)
         aMcumhigh=aMcumhigh+ aMLzhigh(i)
         write (33,341) i, rr*xfact, 
     &      aMcumlow/aMassTotal, aMcumhigh/aMassTotal,aMcum/aMassTotal 
      enddo
 341  format(I4, 4(g14.3))
      end
c-------------------------------------------------------------------- 
c-------------------------------------------------------------------- 
      SUBROUTINE ProfilesSp(xc,yc,zc,T_gas,rmin, rmax,nrbin,fprofname,
     &                      rvirHalo,aMvirHalo)
c-------------------------------------------------------------------- 
c     Make spherical profiles for gas properties and all components 
c     based on Output_Average_profile in ART_Analysis.F
c
c     input in code units.
c-------------------------------------------------------------------- 
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmin, rmax,rlogmin,rlmax
      integer nrbin
      character*256 fprofname 

      integer m
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
      a3 = aexpn**3
      a2 = aexpn**2
      ai = 1.0 / aexpn
      a2i = 1.0 / a2
      gamma1 = gamma - 1.0 
      rho0Cm =1.12d-5*hubble*hubble*Om0*(1.d0 - Y_p) * ai**3 ! o n_H in cm^-3
      rho0C =  rho0 / a3 / 1.E+18     ! in M_sun/pc3 (proper units)
      r0Ch= r0*1000.0 ! in Kpc h-1 (Comoving)
      r0C = r0*1000.0 / hubble  ! in Kpc (Comoving)
      r0P = r0 *1000.0 / hubble * aexpn ! in Kpc (Physical units)
!      deltavir = 119           ! with respect to the matter mean density (code units) at z>>1
!     deltavir = 340. ! same but at z=0
 
      OmegaA= Om0 / ( Om0 + Oml0*a3)
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
            IF (ptemp.lt.T_gas) then
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
      write(12,12) ivirDM, rvirDM*r0P, aMvirDM/hubble, deltavir
c-------------------------------------------------------------------- 
      rvirHalo=rvir*r0P
      aMvirHalo=aMvir/hubble 
c--------------------------------------------------------------------
!     Write outputfile:
      nfpn = index ( fprofname , ' ' ) - 1
      open ( 12 , file = fprofname(1:nfpn) )
      write(12,10) (xc-1.)*r0  , (yc-1.)*r0 , (zc-1.)*r0    
      write(12,11) rmin*r0P, rmax*r0P , nrbin 
      write(12,12) ivir, rvir*r0P, aMvir/hubble, deltavir
      write (12,13) c
      write(12,14) aexpn , 1/aexpn -1. , T_gas
      write(12,15)
      write(12,16)
 10   format ('# (xc,yc,zc)[Mpc h-1 comoving]= ', 3(F9.4,1x))
 11   format ('# rmin = ',g12.6,' rmax =',g12.6,
     &     ' nrbin =',i5)
 12   format ('# ivir = ',i4,' Rvir =',f8.1,
     &     ' [Kpc], Mvir =',g13.6,
     &     ' [Msun] for Dvir(mean)=',f7.2)
 13   format ('NFW fit: concentration= ', f3.0)
 14   format
     &('Distances in Kpc ,densities in cm-3 and M in M_Sun,
     &  for a= ',f6.3, ' and redshift= ', f6.2, 
     &' Baryonic fraction of cold gas', E7.1, ' K ')
 15   format ('# <r>',14x,' r_l',14x,' r_m',14x,' r_r',14x,
     &     ' Rho_g ',3x, 'Rho_St ',8x, 'Rho_dm',10x,
     &     'M_insideR_Gas',3x,  'M_insideR_Stars',3x,
     &     ' M_insideR_DM ',3x,' M_insideR_tot',4x,
     &     'T_g',10x,'S_g',8x,
     &     'N_dm_Pibin',2x,'Baryonic_F  [Uni]' )
 16   format('---------------------')
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
            write(12,17)
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
      close ( 12 )
      write(*,*) 'done writing spherical profiles. exiting...'
      return
      end
c-------------------------------------------------------------------- 
      SUBROUTINE Write_DM_binary(fprofname,Vxc,Vyc, Vzc, xc,yc,zc, rmax,
     & AxX,AxY,AxZ,iopt)
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmax
      real*8 Xph,Yph,Zph,Vxph,Vyph,Vzph
      character*256 fprofname 
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles
!AM distribution:
      real*8 aMassTotal, aMcum, aMcumlow, aMcumhigh
      parameter ( nmaxbins = 1000 )
      parameter (aLzlow=0.1e3)
      real aMtotal(nmaxbins),aMLzlow(nmaxbins), aMLzhigh(nmaxbins)

      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      aL0 = aM0*vfact*xfact
      aL0s=vfact*xfact
      tfact = t0*aexpn**2
      rmax2= rmax**2

!AM distribution:
      drl=rmax/100.
      aMassTotal=0.0
      do i = 1 , nmaxbins
         aMtotal(i)=0.0
         aMLzlow(i)=0.0
         aMLzhigh(i)=0.0
      enddo

      nfn  = index(fprofname, ' ') - 1
      open ( 33, file=fprofname(1:nfn) , form = 'unformatted' )

      Ndm=0
      do is=nsp(1,1),nsp(1,2) ! first DM specie
         Xph =  x(is)-xc
         Yph =  y(is)-yc
         Zph =  z(is)-zc
         Vxph= vx(is)-Vxc
         Vyph=vy(is) -vyc
         Vzph=vz(is) -vzc
         IF (ABS(Xph).lt.rmax.AND.ABS(Yph).lt.rmax.AND.ABS(Zph).lt.rmax)
     &       THEN
            Ndm=Ndm+1
            write(33) is,
     &           Xph*xfact , 
     &           Yph*xfact ,
     &           Zph*xfact ,
     &           vxph*vfact,vyph*vfact,vzph*vfact,
     &           pw(is)*aM0
         ENDIF
      enddo
      write(*,*)"#DM particles inside ",rmax*xfact," Kpc :",Ndm
      write(*,340) pw(nsp(1,1))*aM0*Ndm
 300  format ( 'time=', E20.6, ' yrs , a= ', F10.6)
 301  format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
 302  format (" X Y Z [Kpc] ,Vx  Vy  Vz [ km/s]  M[Msun] age[yr]")
 334  format(I7 , 2x, 8(e14.6,1x))
 335  format(I7 , 2x, 11(e14.6,1x))
 340  format ( ' with total mass equal to ', g10.2, ' Msun')
      close ( 33 ) 
      end
c-------------------------------------------------------------------- 
      SUBROUTINE Write_Stars_binary(fprofname, Vxc, Vyc, Vzc,
     &     xc,yc,zc, rmax,tage,ZDmax,AxX,AxY,AxZ,iopt)
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmax,tage,ZDmax,Vxc, Vyc, Vzc
      real*8 Xph,Yph,Zph, Vxph, Vyph, Vzph
      character*256 fprofname 
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles
      real*8 aMassTotal0
!AM distribution:
      real*8 aMassTotal, aMcum, aMcumlow, aMcumhigh
      parameter ( nmaxbins = 1000 )
      parameter (aLzlow=0.1e3)
      real aMtotal(nmaxbins),aMLzlow(nmaxbins), aMLzhigh(nmaxbins)
      real*8 ageS, ageU, age,tageGyr

      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      aL0 = aM0*vfact*xfact 
      aL0s= vfact*xfact
      tfact = t0*aexpn**2
      rmax2= rmax**2
! Age:
      ageU = age(t)             ! in Gyr
      tageGyr =  tage / 1.e9    ! in Gyr  

!AM distribution:
      drl=rmax/100.
      aMassTotal=0.0
      do i = 1 , nmaxbins
         aMtotal(i)=0.0
         aMLzlow(i)=0.0
         aMLzhigh(i)=0.0
      enddo

      nfn  = index(fprofname, ' ') - 1
      open ( 33, file=fprofname(1:nfn) , form = 'unformatted' )
      write (*,*)  'writing a binary file with the current stellar mass'

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
         
         IF (ABS(Xph).lt.rmax.AND.ABS(Yph).lt.rmax.AND.ABS(Zph).lt.rmax)
     &       THEN
            ageS = ageU - age(dble( tbirth(istar)))
            if (ageS.le.tageGyr) then
               Nstars=Nstars+1
               aMassTotal=aMassTotal +  pw(is)
               aMassTotal0=aMassTotal0 +  pw0(istar)      
               write(33) istar,
     &              Xph*xfact , 
     &              Yph*xfact ,
     &              Zph*xfact ,
     &              Vxph*vfact, Vyph*vfact, Vzph*vfact,
     &              pw(is)*aM0 , 
c     &              pw0(istar)*aM0,
     &              sngl(ageS)
            endif
         ENDIF
      enddo
      write(*,*)"Nstars inside ",rmax*xfact," Kpc :",Nstars
      write(*,*) "Total mass:",  aMassTotal*aM0
      write(*,*) "Total initial mass:", aMassTotal0*aM0
      write(*,340) tage
 300  format ( 'time=', E20.6, ' yrs , a= ', F10.6)
 301  format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
 302  format("X Y Z [Kpc],Vx  Vy Vz [ km/s] M0[Msun] Lx Ly Lz age[Gyr]")
 340  format ( ' younger than ', E6.1, ' yrs')
      close ( 33 ) 
      end
c-------------------------------------------------------------------- 
      SUBROUTINE Write_Stars_binary0(fprofname, Vxc, Vyc, Vzc,
     &     xc,yc,zc, rmax,tage,ZDmax,AxX,AxY,AxZ,iopt)
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmax,tage,ZDmax,Vxc, Vyc, Vzc
      real*8 Xph,Yph,Zph, Vxph, Vyph, Vzph
      character*256 fprofname 
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles
      real*8 aMassTotal0
!AM distribution:
      real*8 aMassTotal, aMcum, aMcumlow, aMcumhigh
      parameter ( nmaxbins = 1000 )
      parameter (aLzlow=0.1e3)
      real aMtotal(nmaxbins),aMLzlow(nmaxbins), aMLzhigh(nmaxbins)
      real*8 ageS, ageU, age,tageGyr

      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      aL0 = aM0*vfact*xfact 
      aL0s= vfact*xfact
      tfact = t0*aexpn**2
      rmax2= rmax**2
! Age:
      ageU = age(t)             ! in Gyr
      tageGyr =  tage / 1.e9    ! in Gyr  

!AM distribution:
      drl=rmax/100.
      aMassTotal=0.0
      do i = 1 , nmaxbins
         aMtotal(i)=0.0
         aMLzlow(i)=0.0
         aMLzhigh(i)=0.0
      enddo

      nfn  = index(fprofname, ' ') - 1
      open ( 33, file=fprofname(1:nfn) , form = 'unformatted' )
      write (*,*)  'writing a binary file with the current stellar mass'

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
         
         IF (ABS(Xph).lt.rmax.AND.ABS(Yph).lt.rmax.AND.ABS(Zph).lt.rmax)
     &       THEN
            ageS = ageU - age(dble( tbirth(istar)))
            if (ageS.le.tageGyr) then
               Nstars=Nstars+1
               aMassTotal=aMassTotal +  pw(is)
               aMassTotal0=aMassTotal0 +  pw0(istar)      
               write(33) istar,
     &              Xph*xfact , 
     &              Yph*xfact ,
     &              Zph*xfact ,
     &              Vxph*vfact, Vyph*vfact, Vzph*vfact,
     &              pw0(istar)*aM0,
     &              sngl(ageS)
            endif
         ENDIF
      enddo
      write(*,*)"Nstars inside ",rmax*xfact," Kpc :",Nstars
      write(*,*) "Total mass:",  aMassTotal*aM0
      write(*,*) "Total initial mass:", aMassTotal0*aM0
      write(*,340) tage
 300  format ( 'time=', E20.6, ' yrs , a= ', F10.6)
 301  format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
 302  format("X Y Z [Kpc],Vx  Vy Vz [ km/s] M0[Msun] Lx Ly Lz age[Gyr]")
 340  format ( ' younger than ', E6.1, ' yrs')
      close ( 33 ) 
      end
c-------------------------------------------------------------------- 
      SUBROUTINE Write_Stars_binaryZ(fprofname, Vxc, Vyc, Vzc,
     &     xc,yc,zc, rmax,tage,ZDmax,AxX,AxY,AxZ,iopt)
      include 'a_tree.h'
      include 'a_control.h'
#    include "a_def.h"

      real xc,yc,zc,rmax,tage,ZDmax,Vxc, Vyc, Vzc
      real*8 Xph,Yph,Zph, Vxph, Vyph, Vzph
      character*256 fprofname 
      real AxX(3), AxY(3), AxZ(3) ! New axes for cylindrical profiles
      real*8 aMassTotal0
!AM distribution:
      real*8 aMassTotal, aMcum, aMcumlow, aMcumhigh
      parameter ( nmaxbins = 1000 )
      parameter (aLzlow=0.1e3)
      real aMtotal(nmaxbins),aMLzlow(nmaxbins), aMLzhigh(nmaxbins)
      real*8 ageS, ageU, age,tageGyr

      xfact = 1000.0*aexpn*r0/hubble 
      vfact = v0 / aexpn 
      aL0 = aM0*vfact*xfact 
      aL0s= vfact*xfact
      tfact = t0*aexpn**2
      rmax2= rmax**2
! Age:
      ageU = age(t)             ! in Gyr
      tageGyr =  tage / 1.e9    ! in Gyr  

!AM distribution:
      drl=rmax/100.
      aMassTotal=0.0
      do i = 1 , nmaxbins
         aMtotal(i)=0.0
         aMLzlow(i)=0.0
         aMLzhigh(i)=0.0
      enddo

      nfn  = index(fprofname, ' ') - 1
      open ( 33, file=fprofname(1:nfn) , form = 'unformatted' )

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
         
         IF (ABS(Xph).lt.rmax.AND.ABS(Yph).lt.rmax.AND.ABS(Zph).lt.rmax)
     &       THEN
            ageS = ageU - age(dble( tbirth(istar)))
            if (ageS.le.tageGyr) then
               Nstars=Nstars+1
               aMassTotal=aMassTotal +  pw(is)
               aMassTotal0=aMassTotal0 +  pw0(istar)      
               write(33) istar,
     &              Xph*xfact , 
     &              Yph*xfact ,
     &              Zph*xfact ,
     &              Vxph*vfact, Vyph*vfact, Vzph*vfact,
     &              pw0(istar)*aM0,
     &              sngl(ageS),
     &              zstII(istar),zstIa(istar)
            endif
         ENDIF
      enddo
      write(*,*)"Nstars inside ",rmax*xfact," Kpc :",Nstars
      write(*,*) "Total mass:",  aMassTotal*aM0
      write(*,*) "Total initial mass:", aMassTotal0*aM0
      write(*,340) tage
 300  format ( 'time=', E20.6, ' yrs , a= ', F10.6)
 301  format ('# (xc,yc,zc)[Mpc h-1]= ', 3(F8.3,1x))
 302  format("X Y Z [Kpc],Vx  Vy Vz [ km/s] M0[Msun] Lx Ly Lz age[Gyr]")
 340  format ( ' younger than ', E6.1, ' yrs')
      close ( 33 ) 
      end


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
      real*8 rs, rsd, rhogl, T_g, Z_met, xHII, xHI, coolrate, cool, heat

!      parameter (QSOBox=28) !300 ) ! 400 ) ! Kpc
      QSOBoxcu = QSOBox/ ( r0 * 1000.0 / hubble * aexpn) ! in code units
      QSOBoxcu2= QSOBoxcu/2.0
      r0Pk= r0*1000.0 / hubble * aexpn ! in Kpc (Physical units)
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

      nfpn = index ( filename , ' ' ) - 1
c
      Ncells=0
c
         open ( 12 , file = filename(1:nfpn), form = 'unformatted' )

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
               write(12) 2.**(-Level)*r0Pk*1000.,
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
                  write(12)   2.**(-Level)*r0Pk*1000.,
     &                 Posxc*r0Pk ,Posyc*r0Pk, Poszc*r0Pk,  
     &                 vx0*v0P,vy0*v0P,vz0*v0P, rho*rho0Cm, ptemp
                  Ncells=Ncells+1

            ENDIF
         endif
      ENDDO
      write (*,*)  'Number of cells= ', Ncells
      CLOSE(12)
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
      real*8 rs, rsd, rhogl, T_g, Z_met, xHII, xHI, coolrate, cool, heat

!      parameter (QSOBox=28) !300 ) ! 400 ) ! Kpc
      QSOBoxcu = QSOBox/ ( r0 * 1000.0 / hubble * aexpn) ! in code units
      QSOBoxcu2= QSOBoxcu/2.0
      r0Pk= r0*1000.0 / hubble * aexpn ! in Kpc (Physical units)
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

      nfpn = index ( filename , ' ' ) - 1
c
      Ncells=0  
c 
         open ( 12 , file = filename(1:nfpn), form = 'unformatted' )

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
                  write(12) 2.**(-Level)*r0Pk*1000.,
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
                  write(12)   2.**(-Level)*r0Pk*1000.,
     &                 Posxc*r0Pk ,Posyc*r0Pk, Poszc*r0Pk,  
     &                 vx0*v0P,vy0*v0P,vz0*v0P, rho*rho0Cm, ptemp,
     &                 ZII, ZIa
                  Ncells=Ncells+1
               ENDIF
            endif
         ENDDO
      write (*,*)  'Number of cells= ', Ncells
      CLOSE(12)
      END
