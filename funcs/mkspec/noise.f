c     PACKAGE: Mockspec
c     PROGRAM: specsynth
c     MODULE:  noise.f
c
c
c     DESCRIPTION
c     routines that insert Gaussian noise into the spectrum based upon
c     the signal to noise ratio given by the user
c
c
c     this file contains:
c     SUBROUTINE addnoise
c     FUNCTION GASDEV
c     FUNCTION ran1
c
c
c..............................................................................
c                                                                             
      SUBROUTINE addnoise(idum)

c     add Gaussian noise to the convolved spectrum, whic up to this time
c     is smooth

c                                                                            
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include           'specsynth.h'
      include           'const.dek'
      integer           pixi,idum
      double precision  xdum,ran1
      double precision  noise,gasdev
      double precision  term,Ic,scale



c     if snr=0 then we are making a noiseless spectrum

      IF (snr.eq.0.0) then
       DO 03 pixi=1,ndata
        sigma(pixi)  = snr
 03    CONTINUE
       RETURN
      END IF

c     use the last random number generating seed to generate a new
c     random seed; this ensures that the noise charactersitics of the
c     spectra have different realizations as a function of pixel number

 01   xdum = ran1(idum)
      If (xdum.lt.0.1) GOTO 01
      idum = int(-1*xdum*getpid())
c      idum = int(-50.0*xdum)

c     add the noise and make the sigma spectrum; uses the formula from
c     my thesis in which the SNR is the critical quantity.  the SNR is
c     used tom estimate the continuum level, Ic, which provides the
c     scale of the deviations using a Poisson model for large number
c     statistics (basically a Gaussian distribution)

       term = 1.0d0 + 4.0d0*(rn/snr)**2
       Ic   = 0.5d0*snr**2 * (1.0d0+sqrt(term))
       DO 05 pixi=1,ndata
        scale        = sqrt(wrkflx(pixi)*Ic + rn**2)/Ic
        noise        = scale*gasdev(idum)
        wrkflx(pixi) = wrkflx(pixi) + noise
        sigma(pixi)  = scale
 05   CONTINUE

      RETURN

      END

c
c.......................................................................
c

      DOUBLE PRECISION FUNCTION GASDEV(IDUM)

c     Numerical Recipe code: uses function RAN1 to generate a random
c     Gaussian deviate deviate (weighted by the cumulative distribution
c     function of a Gaussian function

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE ISET,GSET
      DATA ISET/0/

      IF (ISET.EQ.0) THEN
1       V1=2.*RAN1(IDUM)-1.
        V2=2.*RAN1(IDUM)-1.
        R=V1**2+V2**2
        IF(R.GE.1.) GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        ISET=1
      ELSE
        GASDEV=GSET
        ISET=0
      ENDIF

      RETURN
      END

c
c.......................................................................
c

      DOUBLE PRECISION FUNCTION ran1(idum)

c     Numerical Recipe code: generates a random deviate (real number
c     between 0 and 1) given an integer seed "idum" for the kernal

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
C      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     $     NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=3.e-16,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/

      if (idum.le.0.or.iy.eq.0) then
         idum=max(-idum,1)
         do j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
         enddo
         iy=iv(1)
      endif

      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)

      return
      END
