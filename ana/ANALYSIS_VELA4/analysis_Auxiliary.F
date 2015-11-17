C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
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
