c     PACKAGE: Mockspec/anal
c     PROGRAM: sysanal
c     MODULE:  recipes.f
c
c     DESCRIPTION
c     this module contains modified versions of Numerical Recipe
c     Routines
c
c
c     this file contains:
c     subroutine   spline
c     subroutine   splint
c     subroutine   qromb
c     subroutine   trapzd
c     subroutine   polint
c     subroutine   indexx
c     dbl function value2, which is not from the NR library
c
c
c..........................................................................
c
      subroutine spline(x,y,n,y2)   
c
c     given arrays x and y of length n containing a tabulated function
c     i.e y=f(x), with the x monotonically increasing and given values
c     for yp1 and ypn, the first derivative at the points 1 and n
c     respectively, this routine returns the array y2 of length n which
c     contains the second derivatives of the at the tabulated points x.
c     if yp1 and/or ypn are set larger than 1e30, the routine sets the
c     boundary condtion for a natural spline (one with zero second
c     derivative) at the boundary.
c
c     this routine is only called once for any given x and y.
c 
c..........................................................................

      implicit none
      save
      integer            n,i,k,nmax 
      parameter          (nmax=10000) 
      double precision   x(nmax),y(nmax),y2(nmax),yp1,ypn,u(nmax),
     &                   sig,p,qn,un
      parameter          (yp1=1.0e33, ypn=1.0e33)


c     the lower boundary condition is set to be either natural

      if (yp1 .gt. .99E30) then 
        y2(1)=0.0   
        u(1)=0.0

c     or it has a specified first derivative

      else  
        y2(1)=-0.5  
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      end if

c     this is the decomposition loop of the tridiagonal algorithm. y2
c     and u are used for temporary storage of the decomposition factors.

      do 11 i=2,n-1 
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))   
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) 
     +      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p   
11    continue  
 
c     the upper boundary condition is set to be either natural

      if (ypn .gt. .99E30) then 
        qn=0.0  
        un=0.0  
 
c     or it has a specified first derivative

      else  
        qn=0.5  
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))) 
      end if

c     this is the backsubstitution loop of the tridiagonal algorithm

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)  
      do 12 k=n-1,1,-1  
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue  
 
c     return

      return
      end   

c
c..........................................................................
c 
      subroutine splint(xa,ya,y2a,n,x,y)

c
c     given the arrays xa and ya of length n, which tablulate a
c     monotonic function, and given the array y2a, which is the output
c     of spline (above), and given a value of x this routine returns a
c     cubic spline interpolated value of y.
c 
c..........................................................................


      implicit none
      save
      integer             n,klo,khi,k   
      double precision    xa(n),ya(n),y2a(n),x,y,h,a,b  


c     find the right place in the table by bisection. this is optimal if
c     the sequential calls to this routine are at random values of x.
c     if the sequential calls are in order and closely spaced, one might
c     store the values of klo and khi and test if they remain
c     appropriate on next call

      klo=1 
      khi=n 
1     if (khi-klo .gt. 1) then  
       k=(khi+klo)/2
       if (xa(k) .gt. x) then   
        khi=k   
       else 
        klo=k   
       end if   
       goto 1   
      end if    

c     klo and khi now bracket the input value of x the xa's must be
c     distinct

      h=xa(khi)-xa(klo) 
      if (h .eq. 0.0d0) then
        write(6,99) khi,xa(khi),klo,xa(klo)
   99   format(1x,'  khi=',i4,'  xa(khi)=',1pe13.5,'  klo=',
     1         i4,'  xa(klo)=',1pe13.5)
        stop 'bad xa input in routine splint' 
      end if

c     evaluate the cubic spline

      a=(xa(khi)-x)/h   
      b=(x-xa(klo))/h   
      y=a*ya(klo)+b*ya(khi)+
     +      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0

c     return

      return    
      end   

c
c..........................................................................
c 
      subroutine qromb(func,a,b,eps,ss) 

c 
c     returns as ss the integral of the function func from a to b with
c     fractional accuracy eps. integration by romberg's method of order
c     2k where e.g k=2 is simpson's rule.
c 
c     jmax limits the total number of steps; k is the the number of
c     points used in the extrapolation; arrays s and h store the
c     trapazoidal approximations and their relative step sizes.
c
c..........................................................................

      implicit none
      save
      external          func
      integer           j,jmax,jmaxp,k,km   
      parameter         (jmax=20, jmaxp=jmax+1, k=5, km=k-1)
      double precision  a,b,ss,s(jmaxp),h(jmaxp),eps,dss,zero 



      zero = 0.0d0
      h(1) = 1.0

      do 11 j=1,jmax
C       WRITE(6,*) 'call trapzd'
       call trapzd(func,a,b,s(j),j) 
       if (j .ge. k) then   
C       WRITE(6,*) 'call polint'
        call polint(h(j-km),s(j-km),k,zero,ss,dss)   
        if (abs(dss) .le. eps*abs(ss)) return   
       end if   
       s(j+1) = s(j)
       h(j+1) = 0.25 * h(j) 
11    continue  

      write(6,*) ' after ',jmax,' iterations '
      write(6,*) ' of trying to integrate between ',a,' and ',b
      write(6,*) ' and fractional accuracy ',eps
      write(6,*) ' the integral is ',ss
      write(6,*) ' and error estimate ',dss
      write(6,*) ' so that abs(dss) ',abs(dss),
     +           ' > eps*abs(ss)',eps*abs(ss)
      stop       'too many steps in qromb' 

c     return

      end   
 
c 
c..........................................................................
c 
      subroutine trapzd(func,a,b,s,n)   
c 
c     this routine computes the n'th stage of refinement of an extended
c     trapazoidal rule. func is input as the name of a function to be
c     integrated between limits a and b. when n=1 the routine returns as
c     s the crudest estimate of the integral of func(x)dx from a to b.
c     subsequent calls with n=2,3... will improve the accuracy of s by
c     adding 2**(n-2) additional interior points. s should not be
c     modified between sequential calls.
c
c     local it is the number of points to be added on the next call
c     local del is the step size.
c
c..........................................................................

      implicit none
      save
      external          func
      integer           n,it,j  
      double precision  func,a,b,s,del,x,sum,tnm
 
c     go

      if (n.eq.1) then  
       s = 0.5 * (b-a) * ( func(a) + func(b) )  
       it = 1   
      else  
       tnm = it 
       del = (b-a)/tnm  
       x = a + (0.5 *del)   
       sum = 0.0
       do 11 j=1,it 
        sum = sum + func(x) 
        x = x + del 
11     continue 
       s = 0.5 * (s + (b-a)*sum/tnm)
       it = 2 * it  
      end if

      return    
      end   

c
c..........................................................................
c

      subroutine polint(xa,ya,n,x,y,dy)

c
c     given arrays xa and ya of length n and a value x, this routine
c     returns a value y and an error estimate dy. if p(x) is the
c     polynomial of degree n-1 such that ya = p(xa) ya then the returned
c     value is y = p(x)
c
c..........................................................................


      implicit none
      save
      integer          n,nmax,ns,i,m
      parameter        (nmax=30)
      double precision xa(1),ya(1),x,y,dy,
     +                 c(nmax),d(nmax),dif,dift,ho,hp,w,den

c     find the index ns of the closest table entry and initialize the c
c     and d tables

      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue

c     first guess for y

      y=ya(ns)

c     for each column of the table, loop over the c's and d's and update
c     them

      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) stop ' 2 xa entires are the same in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue

c     after each column is completed, decide which correction c or d, to
c     add to the accumulating value of y, that is, which path to take in
c     the table by forking up or down. ns is updated as we go to keep
c     track of where we are. the last dy added is the error indicator.

        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue

      return
      end

c..........................................................................
c 
      double precision function zbrent(func,x1,x2,tol)  

c 
c     using brent's method, find the root of a function func known to
c     lie between x1 and x2. the root is returned as zbrent and refined
c     until its accuracy is tol
c..........................................................................
c

      implicit none
      save
      integer           itmax,iter  
      parameter         (itmax=100)  
      double precision  func,x1,x2,tol,a,b,c,d,e,fa,
     1                  fb,fc,xm,tol1,p,q,r,s,eps   
      parameter         (eps=3.0e-8)  

c     initialize

      a  = x1
      b  = x2
      fa = func(a)  
      fb = func(b)  
      if (fb*fa .gt. 0.0d0) then
       write(6,99) a,fa,b,fb
 99    format(1x,'  a=',1pe13.5,'  f(a)=',1pe13.5,
     1        '  b=',1pe13.5,'  f(b)=',1pe13.5) 
       stop 'root not bracketed in routine zbrent'   
      end if
      fc = fb   

      do 11 iter =1,itmax   

c     rename a,b,c and adjusting bound interval d

       if (fb*fc .gt. 0.0) then 
        c  = a   
        fc = fa 
        d  = b-a 
        e  = d   
       end if   
       if (abs(fc) .lt. abs(fb)) then   
        a  = b   
        b  = c   
        c  = a   
        fa = fb 
        fb = fc 
        fc = fa 
       end if   
       tol1 = 2.0d0 * eps * abs(b) + 0.5d0 * tol
       xm   = 0.5d0 * (c-b) 

c     convergence check

       if (abs(xm) .le. tol1 .or. fb .eq. 0.0) then 
        zbrent = b  
        return  
       end if   

c     attempt quadratic interpolation

       if (abs(e) .ge. tol1 .and. abs(fa) .gt. abs(fb)) then
        s = fb/fa   
        if (a .eq. c) then  
         p = 2.0d0 * xm * s   
         q = 1.0d0 - s 
        else
         q = fa/fc  
         r = fb/fc  
         p = s * (2.0d0 * xm * q *(q-r) - (b-a)*(r - 1.0d0))  
         q = (q - 1.0d0) * (r - 1.0d0) * (s - 1.0d0)
        end if  

c     check if in bounds

        if (p.gt. 0.0) q = -q   
        p = abs(p)  

c     accept interpolation

        if (2.0d0*p .lt. min(3.0d0*xm*q - abs(tol1*q),abs(e*q))) then   
         e = d  
         d = p/q

c     or bisection

        else
         d = xm 
         e = d  
        end if  

c     bounds decreasing to slowly use bisection

       else 
        d = xm  
        e = d   
       end if   

c     move best guess to a

       a  = b
       fa = fb  
       if (abs(d) .gt. tol1) then   
        b = b + d   
       else 
        b = b + sign(tol1,xm)   
       end if   
       fb = func(b) 
11    continue  

      stop 'too many iterations in routine zbrent'  

c     done

      end   


c..........................................................................
c

      subroutine indexx(n,arrin,indx)   

c
c     indexes an array arrin of length n. outputs the array indx such
c     that arrin(index(j)) (j=1,...n) is in ascending order. the input
c     array arrin is not changed.
c
c..........................................................................


      implicit none
      save
      integer          n,indx(n),indxt,j,l,ir,i
      double precision arrin(n),q

c
c     iniitalize the index array with consecutive integers

      do 11 j=1,n   
        indx(j)=j   
11    continue  

c     from here on out its heapsort, but with indirect indexing through
c     indx in all references to arrin. compare to the heapsort routine
c     above.

      l=n/2+1   
      ir=n  

10    continue  

        if(l.gt.1)then  
          l=l-1 
          indxt=indx(l) 
          q=arrin(indxt)
        else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)  
          ir=ir-1   
          if(ir.eq.1)then   
            indx(1)=indxt   
            return  
          endif 
        endif   

        i=l 
        j=l+l   

20      if(j.le.ir)then 
          if(j.lt.ir)then   
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1 
          endif 
          if(q.lt.arrin(indx(j)))then   
            indx(i)=indx(j) 
            i=j 
            j=j+j   
          else  
            j=ir+1  
          endif 
        go to 20
        endif   
        indx(i)=indxt   

      go to 10  

      end   


c..........................................................................
c

      double precision function value2(string,err)

c
c     this routine takes the character string number and converts it to
c     a real. if trouble is encountered during the conversion, this
c     routine returns the logical err as .true.
c
c     this is not a Numerical Reicpe, it is home grown
c
c..........................................................................


      implicit none
      save
      logical          pflag,err
      character*(*)    string
      character*1      plus,minus,decmal,blank,se,sd,se1,sd1
      integer          noblnk,long,ipoint,power,psign,iten,j,z,i
      double precision x,sign,factor,rten,temp
      parameter        (plus = '+'  , minus = '-' , decmal = '.'   ,
     1                  blank = ' ' , se = 'e'    , sd = 'd'       ,
     2                  se1 = 'E'   , sd1 = 'D'   , rten =  10.0,
     3                  iten = 10                                   )


c     initialize

      err    =  .false.
      x      =  0.0d0
      sign   =  1.0d0
      factor =  rten
      pflag  =  .false.
      noblnk =  0
      power  =  0
      psign  =  1
      long   =  len(string)

c     remove leading blanks and the sign

      do 10 z = 1,7
       noblnk = noblnk + 1
       if ( string(noblnk:noblnk) .eq. blank) then
        if (noblnk .gt. 6 ) goto 1000
       else
        if (string(noblnk:noblnk) .eq. plus) then
         noblnk = noblnk + 1
        else if (string(noblnk:noblnk) .eq. minus) then
         noblnk = noblnk + 1
         sign =  -1.0d0
        end if
        go to 100
       end if
10    continue

c     main number conversion loop

100   do 200 i = noblnk,long
       ipoint = i + 1

c     if blank character then we are done

       if ( string(i:i) .eq. blank ) then
        x = x * sign
        value2 = x 
        return
c     if it is an exponent process it

       else if (string(i:i).eq.se  .or. string(i:i).eq.sd .or.
     1          string(i:i).eq.se1 .or. string(i:i).eq.sd1   ) then
        if (x .eq. 0.0d0 .and. ipoint.eq.2)     x = 1.0d0
        if (sign .eq. -1.0d0 .and. ipoint.eq.3) x = 1.0d0
        if (string(ipoint:ipoint) .eq. plus) ipoint = ipoint + 1
        if (string(ipoint:ipoint) .eq. minus) then
         ipoint = ipoint + 1
         psign = -1
        end if
        do 150 z = ipoint,long
         if (string(z:z) .eq. blank)  then
          x = sign * x * rten**(power*psign)
          value2 = x
          return
         else
          j = ichar(string(z:z)) - 48
          if ( (j.lt.0) .or. (j.gt.9) ) go to 1000
          power= (power * iten)  + j
         end if
150     continue

c     if it is a number process it

       else if (string(i:i) .ne. decmal) then
        j = ichar(string(i:i)) - 48
        if ( (j.lt.0) .or. (j.gt.9) ) go to 1000
        if (.not.(pflag) ) then
         x = (x*rten) + j
        else
         temp   = j
         x      = x + (temp/factor)
         factor = factor * rten
         go to 200
        end if

c     must be a decimal point

       else
        if (pflag) go to 1000
        pflag = .true.
       end if
200   continue

c     error trap

1000  err = .true.

c     return

      return
      end

c     eof
