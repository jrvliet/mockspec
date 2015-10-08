c     PACKAGE: Mockspec
c     PROGRAM: specsynth
c     MODULE:  instrument.f
c
c
c     DESCRIPTION
c     the routines herein set up the instrumental spread function for
c     convolving the original model spectrum
c
c     this file contains:
c     subroutine instrument
c     dble function phi
c
c
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
c
      subroutine          instrument(m,wcen,flag)

c     given the instrumental profile sigma in velocity units; the ISF
c     profiles are loaded in wrap-around order for the convolution; zero
c     spatial information is the first index

c     flag = 0 ; call for setting up convolution
c     flag = 1 ; =0 + communicates instrumental parameters

c
c.......................................................................
c

      include             'specsynth.h'
      include             'const.dek'

      integer             i,m,iplus,iminus,np2,flag
      parameter           (np2 = 17)
      double precision    xdv,norm,phi,wcen
      integer             pwrsof2(np2)

      data pwrsof2 /16,32,64,128,256,512,1024,2048,4096,8192,
     &              16384,32768,65536,131072,262144,524288,
     &              1048576/




c     compute the instrumental resolution in velocity units [km/s] sigma
c     in km/s given by FWHM/2.35

      if (R_fac.eq.0.0d0) then
       profile = slit
      else
       profile = slit*ckms/(2.35*R_fac)
      end if

c     dv is the velocity sampling of the pixels km/s/pixel...  the
c     number of pixels per resolution element = profile/dv

      dv = dwave/wcen * ckms
      pixpres = 2.35 * profile/(dv*slit)
      hdv = dv/resfac
      nresponse  = 2*int(conwindo*profile/hdv) + 1

c      write(*,*) "Instrument"
c      write(*,*) profile, dwave, wcen, dv, pixpres, hdv, nresponse
c      write(*,*) dv, resfac, hdv
c     now stuff the response function in wrap around order

      response(1) = phi(0.0d0,profile)
      norm        = response(1)
      do 11 i=1,int(nresponse/2)
       xdv              = real(i)*hdv 
       iplus            = i + 1
       iminus           = nresponse - (i-1)
       response(iplus)  = phi(xdv,profile)
       response(iminus) = response(iplus)
       norm = norm + response(iplus) + response(iminus)
 11   continue
  
c     for the convolution integral, the response function needs to be
c     normalized or flux is not conserved...  unit width is used because
c     data are discretized by pixel indices in the convolution

      do 13 i=1,nresponse
       response(i) = response(i)/norm
 13   continue

      open(unit=82, file="fresponse.dat")
      do 83 i=1,nresponse
       write(82,*) response(i)
 83   continue
      close(83)



c     compute the length of the convolution functions..  NCONDAT and
c     NFFT; do a table search on the powers of 2; this ensures the
c     smallest padding of the CONVDATA array in routine convolve

      ncondat = int(resfac)*(m-1) - 1
      nfft   = ncondat + (nresponse-1)/2 + 1 
      do 15 i=1,np2
       if (nfft.le.pwrsof2(i)) then
        nfft = pwrsof2(i)
        GOTO 17
       end if
 15   continue

c     this is seen only if we exceed the powers of 2 table

      write(STDOUT,*) ' NRESPONSE = ',nresponse
      write(STDOUT,*) ' M         = ',m
      write(STDOUT,*) ' NCONDAT   = ',ncondat
      write(STDOUT,*) ' NFFT      = ',nfft
      write(STDOUT,*) ' NFFT MAX  = ',maxcon
      stop ' ERROR:(instrument): NFFT not defined? too big?'

 17   if (nfft.gt.maxcon) then
       stop ' ERROR:(instrument): nfft > maxcon'
      end if

c     debug, set flag=1 in calling routine

      if (flag.eq.1) then 
       if (R_fac.eq.0.0d0) then
        write(STDOUT,'(a,f8.3)') 
     @  ' Instrumental Sigma     [km/s] =',slit
       else
        write(STDOUT,'(a,f8.0)') 
     @  ' Spectrograph R     [lam/dlam] =',R_fac
        write(STDOUT,'(a,f8.3)') 
     @  ' Effective Slit Width [arcsec] =',slit
        write(STDOUT,'(a,f8.3)') 
     @  ' Instrumental Sigma     [km/s] =',profile
       end if
       write(STDOUT,'(a,f8.3)') 
     @  ' Pixel Resolution       [km/s] =',dv
       write(STDOUT,'(a,f8.3)') 
     @ ' Pix / Res Elem         [FWHM] =',pixpres
       write(STDOUT,'(a,f8.3)') 
     @ ' Response Widths       [sigma] =',conwindo
       write(STDOUT,'(a,f8.3)') 
     @ ' Convolution Res       [1/pix] =',resfac
       write(STDOUT,'(a,i5)')   
     @ ' Response Length         [pix] =',nresponse
       write(STDOUT,'(a,i5)')   
     @ ' Convolution Pixels      [pix] =',nfft
      end if

  
      return

      end


c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  
  
      double precision function phi(dvw,width)
  
c    given the sigma in velocity units this routine computes the 
c    relative Gaussian value of the instrumental profile.  
c    called by routine instrument iteratively

c
c.......................................................................
  
      include             'specsynth.h'
      double precision     dvw,z,width



c     dv is value at which the instrumental response is to be evaluated

      z   = dvw/width
      phi = dexp(-(z**2)/2.0d0)

      return
      end
