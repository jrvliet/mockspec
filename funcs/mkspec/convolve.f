c     PACKAGE: Mockspec
c     PROGRAM: specsynth
c     MODULE:  convolve.f
c
c
c     DESCRIPTION
c     herein, the model spectrum is convolved with the instrumental
c     spread function to mkae the resulting model spectrum (prior to
c     adding noise)
c
c
c..............................................................................
c  

      SUBROUTINE          convolve
 
c     this routine sets up the call to the Num Recipes convolver
c     convolution resolution is higher than data resolution so that it
c     data are "smoother"

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include             'specsynth.h'
      integer             i,k,nmax2,ndat,pixi,mfft
      parameter           (nmax2 = maxcon*2)
      double precision    xa(maxvec),ya(maxvec),x,y,y2a(maxvec)
      double precision    yp1,ypn
      double precision    ans(nmax2+2),phiov(maxcon)

       ndat    = ndata
       mfft    = nfft
       yp1     = 1.0d33
       ypn     = yp1
      


c     stuff the "smoothing" data array, one only need keep the index as
c     the abscissa; then obtain the spline coefficients

      do 27 pixi=1,ndata
       xa(pixi)  = real(pixi)
       ya(pixi)  = wrkflx(pixi)
       y2a(pixi) = 0.0d0
 27   continue
      call spline(xa,ya,ndat,yp1,ypn,y2a)
  
c     stuff the end points and the interpolation points

      convdata(1)       = ya(1)
      convdata(ncondat) = ya(ndat)
      do 29 i=2,ncondat-1
       x = 1.0d0 + real(i-1)/resfac
       call splint(xa,ya,y2a,ndat,x,y)
       convdata(i) = y
 29   continue

c     pad the data array with unity, the continuum flux level this
c     assumes that no features are on the edges of the data array

      do 25 i=ncondat+1,nfft
       convdata(i) = 1.0d0
 25   continue
      
      write(*,*) "nfft   ncondat"
      write(*,*) nfft, ncondat
c     the FFT's do violence to the response function- so re-load it into
c     phi(v)

      do 26 i=1,nresponse
        phiov(i) = response(i)
 26   continue

c     we are "go" for the FFT's

      CALL convlv(convdata,mfft,phiov,nresponse,+1,ans)

c     now, restuff the wrkflx array and bail; pick off every resfac
c     element and stuff into idx element of wrkflx

        do 33 pixi=1,ndata
         k = 1 + (pixi-1)*int(resfac)
         wrkflx(pixi) = abs(ans(k))
 33     continue

      return
      end
