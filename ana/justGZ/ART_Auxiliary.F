c     ===================================================================
c
c      ART Version 4: ART_Auxiliary.f - auxiliary routines & functions
c
c     ===================================================================
c   
c     Routines from the Numerical Recipe : 
c      DTDA
c      INTEGRATE
c      RUNGE5VAR
c      RUNGE
c      ZBRENT
c      RANDD
c
c     ------------------------------------------
      double precision function DTDA(a, fp, np )
c     ------------------------------------------
c
c     Find dt/da given a, from the Friedmann Equation.  This is exact for
c     any isotropic-metric cosmology consistent with General Relativity.
c     Here, "t" is understood to be in units of the inverse Hubble constant
c     (i.e. "t" = H0*t).
c
c     Definitions for parameters are as in Peebles 1993, eqn (5.53).
c
c     a - expansion factor
c     fp(1) - Om0 - matter density wrt critical at z=0
c     fp(2) - Oml0 
c
      implicit none
      
      integer np 
      double precision a, fp(np), omegam, omegal

      double precision temp

      omegam = fp(1) 
      omegal = fp(2)

      temp = omegam + omegal*a**3 + (1.d0-omegam-omegal)*a

      if (temp.gt.0.d0) then
        DTDA = sqrt(a/temp)
      else
        DTDA = 0.d0
      end if

      return
      end

c     --------------------------------------------------------------
      double precision function INTEGRATE(FUNC,fp,np,a,b,dxinit,eps)
c     --------------------------------------------------------------
c
c     Quadrature using fifth order Runge-Kutta with adaptive step size.
c     Based on Press et al, Numerical Recipes in C, 2nd ed, pp 719-722.
c
c     Runge-Kutta driver with adaptive stepsize control.  Integrate starting
c     value y from a to b with accuracy eps, storing intermediate results in
c     global variables.  dxinit should be set as a guessed first stepsize.
c
c     Pass a vector of parameters of length np to FUNC in fp(np).
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC

      implicit none
      integer np 
      double precision a, b, eps, dxinit, FUNC, fp(np)
      external FUNC

      integer maxsteps
      parameter(maxsteps=100000000)

      double precision x, dx, dxnext, y, dydx, yscale
      integer  Nstep

      x     = a
      dx    = dxinit
      y     = 0.d0
      Nstep = 0

      do while ((x-b)*(b-a).lt.0.d0.and.Nstep.lt.maxsteps)
        Nstep = Nstep + 1
        dydx = FUNC(x,fp,np)
c
c       yscale is the scaling used to monitor accuracy.  This general-purpose
c       choice can be modified if need be.
c
        yscale = max(abs(y) + abs(dx*dydx), 1.d-20)
        if ((x+dx-b)*(x+dx-a).gt.0.d0)  ! If stepsize overshoots, decrease it.
     1    dx = b - x

        call RUNGE5VAR(y,dydx,x,dx,eps,yscale,dxnext,FUNC,fp,np)

        dx = dxnext
      end do

      if (Nstep.ge.maxsteps)
     1  write (*,*) 'WARNING: failed to converge in INTEGRATE.'

      INTEGRATE = y

      return
      end
c
c     -------------------------------------------------------------
      SUBROUTINE RUNGE5VAR(y,dydx,x,htry,eps,yscale,hnext,DERIVS,
     1                     fp,np)
c     -------------------------------------------------------------
c
c     Fifth-order Runge-Kutta step with monitoring of local truncation error
c     to ensure accuracy and adjust stepsize.  Input are the dependent
c     variable y and its derivative dydx at the starting value of the
c     independent variable x.  Also input are the stepsize to be attempted
c     htry, the required accuracy eps, and the value yscale, against which the
c     error is scaled.  On output, y and x are replaced by their new values.
c     hdid is the stepsize that was actually accomplished, and hnext is the
c     estimated next stepsize.  DERIVS is the user-supplied routine that
c     computes right-hand-side derivatives.  The argument fp is a vector 
c     of parameters (np parameters) to be passed to DERIVS 
c     (NOT integrated over).
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC
c
      implicit none
      integer np 
      double precision eps,hnext,htry,x,dydx,y,yscale,DERIVS,fp(np)
      external DERIVS

      double precision errmax,h,hold,htemp,xnew,yerr,ytemp

      double precision safety,pgrow,pshrink,errcon
      parameter (safety  =  0.9d0)
      parameter (pgrow   = -0.2d0)
      parameter (pshrink = -0.25d0)
      parameter (errcon  =  1.89d-4)

      h = htry                         ! Set stepsize to initial accuracy.
      errmax = 10.d0
      do while (errmax.gt.1.d0)
        call RUNGE(y,dydx,x,h,ytemp,yerr,DERIVS,fp,np)

        errmax = abs(yerr/yscale)/eps   ! Scale relative to required accuracy.
        if (errmax.gt.1.d0) then        ! Truncation error too large; reduce h
          htemp = safety*h*(errmax**pshrink)
          hold = h
          h = sign(max(abs(htemp),0.1d0*abs(h)),h)  ! No more than factor of 10
          xnew = x + h
          if (xnew.eq.x) then
            write (*,*) 'WARNING: ',
     1                  'Stepsize underflow in RUNGE5VAR().'
            h = hold
            errmax = 0.d0
          end if
        end if
      end do
c
c     Step succeeded.  Compute estimated size of next step.
c
      if (errmax.gt.errcon) then
        hnext = safety*h*(errmax**pgrow)
      else
        hnext = 5.d0 * h                ! No more than factor of 5 increase.
      end if
      x = x + h

      y = ytemp

      return
      end
c
c     ---------------------------------------------------
      SUBROUTINE RUNGE(y,dydx,x,h,yout,yerr,DERIVS,fp,np)
c     ---------------------------------------------------
c
c     Given values for a variable y and its derivative dydx known at x, use
c     the fifth-order Cash-Karp Runge-Kutta method to advance the solution
c     over an interval h and return the incremented variables as yout.  Also
c     return an estimate of the local truncation error in yout using the
c     embedded fourth order method.  The user supplies the routine
c     DERIVS(x,y,dydx), which returns derivatives dydx at x.
c
c     Copyright (c) 1997 Michael A. K. Gross.  You may use this program for
c     personal, educational or research purposes.  Commercial purposes require
c     special arrangements. If you publish a paper that depends upon this code,
c     please cite it appropriately.
c
c     Questions and/or comments may be sent to gross@fozzie.gsfc.nasa.gov.
c
c     slight modification by A.Kravtsov to allow to pass a vector 
c     of parameters to FUNC
c
      implicit none

      integer np
      double precision h,x,dydx,y,yerr,yout,DERIVS,fp(np)

      external DERIVS

      double precision ak3, ak4, ak5 ,ak6

      double precision a2,a3,a4,a5,a6
      double precision c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
      parameter(a2  =    0.2d0)
      parameter(a3  =    0.3d0)
      parameter(a4  =    0.6d0)
      parameter(a5  =    1.d0)
      parameter(a6  =    0.875d0)
      parameter(c1  =   37.d0/378.d0)
      parameter(c3  =  250.d0/621.d0)
      parameter(c4  =  125.d0/594.d0)
      parameter(c6  =  512.d0/1771.d0)
      parameter(dc1 = c1 -  2825.d0/27648.d0)
      parameter(dc3 = c3 - 18575.d0/48384.d0)
      parameter(dc4 = c4 - 13525.d0/55296.d0)
      parameter(dc5 = -277.d0/14336.d0)
      parameter(dc6 = c6 -     0.25d0)

      ak3 = DERIVS(x+a3*h,fp,np)
      ak4 = DERIVS(x+a4*h,fp,np)
      ak5 = DERIVS(x+a5*h,fp,np)
      ak6 = DERIVS(x+a6*h,fp,np)
c
c     Estimate the fifth order value.
c
      yout = y + h*(c1*dydx + c3*ak3 + c4*ak4  + c6*ak6)
c
c     Estimate error as difference between fourth and fifth order
c
      yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)

      return
      end
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
        write(*,*) 'PAUSE: Code breaks.'
        read(*,*)
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
      write(*,*) 'PAUSE: Code breaks. ',
     & 'zbrent exceeding maximum iterations'
      read(*,*)
      zbrent=b
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
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
