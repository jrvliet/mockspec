
      DOUBLE PRECISION FUNCTION phfit2(nz,ne,is,e)

CCC Verner's original comments (though in press paper was updated)
*********************************************************************************
*** Version 2. March 25, 1996.
*** Written by D. A. Verner, verner@pa.uky.edu
*** Inner-shell ionization energies of some low-ionized species are slightly
*** improved to fit smoothly the experimental inner-shell ionization energies 
*** of neutral atoms.
*********************************************************************************
*** This subroutine calculates partial photoionization cross sections
*** for all ionization stages of all atoms from H to Zn (Z=30) by use of
*** the following fit parameters:
***
*** Outer shells of the Opacity Project (OP) elements:
*** Verner, Ferland, Korista, Yakovlev, 1996, ApJ, 465, 487
***
*** Inner shells of all elements, and outer shells of the non-OP elements:
*** Verner and Yakovlev, 1995, A&AS, 109, 125
***
*** Input parameters:  nz - atomic number from 1 to 30 (integer) 
***                    ne - number of electrons from 1 to iz (integer)
***                    is - shell number (integer)
***                    e - photon energy, eV 
***
*** Output parameter:  s - photoionization cross section, Mb
***
*** Shell numbers:     1=1s, 2=2s, 3=2p, 4=3s, 5=3p, 6=3d, 7=4s. 
***
*** If a species in the ground state has no electrons on the given shell,
*** the subroutine returns s=0.
*********************************************************************************

c     all comments below this point are by Chris Churchill; they make a
c     lot more sense if you have the papers in front of you!
c
c     note that Verner's original code is a subroutine, I have modified
c     this to be a function, and thus have removed the call parameter
c     "s", since "phfit2" will be set to "s" prior to return
c
c     **** NOTE ***
c     Verner's code bombs for (iz,ne)=(19,19),(20,19),(20,20) because
c     the 3d shells are empty; the 4s shell fills before the 3d shell
c     (Hund's rule); I thus modified the code to skip s=6 for these ions
c
c     the fitting formula is
c
c     s(e) = s0*F(y)  [Mb = 1.0e-18 cm^2]
c
c     where
c
c     F(y) = [(x-1)^2 + yW^2] * (y^q) * [1+sqrt(y/yA)]^-P
c
c     where
c
c     x = e/e0 - y0
c     y = sqrt(x^2 + y1^2)
c  *  q = 0.5*P - 5.5
c
c     and where the fitting parameters are:
c     s0(nz,ne) = normalization constant, s(e) amplitude scale
c     e0(nz,ne) = energy scale
c     yW(nz,ne) = provides s(e) dip energy width (gamma=yW*e0)
c     yA(nz,ne) = provides upper energy boundry of powerlaw (eA=yA*e0)
c  ** P(nz,ne)  = powerlaw energy scaling of xsec (energy range e2<e<eA)
c     y0(nz,ne) = provides xsec maximum (max=y0*e0)
c     y1(nz,ne) = provides xsec energy width (gamma=y1*e0) 
c
c  -  *  for subshell L(is), q = 0.5*P -5.5 - L(is) , to obtain correct
c        asymptotic behavior
c  -  ** e2=e1+e0, is the position of s(e) energy dip (minimum), where 
c        e1=y0*e0
c
c     note that the powerlaw, q, is valid between the energy dip, e2,
c     and the upper boundary energy, eA
c
c     as stated above
c     nz = species atomic number ( = k )
c     ne = number of electrons for ion ( = k-j+1, where j = ion stage )
c     is = shell number electron is in
c
c     INNER SHELLS: Verner uses table parameters PH1(I,nz,ne,is) with
c     eX(I=1) e0(I=2) s0(I=3) yA(I=4) P(I=5) yW(I=6) 
c
c     if is=ntot(I) eX=eth     (threshold outer e-) 
c     if is=1       eX=Emax    (threshold most tightly bound e-)
c     otherwise     eX=eth(is) (threshold e- in shell=is)
c     note that in the parameter set y0=0 for all entries
c
c     OUTER SHELL: Verner uses table parameters PH2(I,nz,ne) with
c     e0(I=1) s0(I=2) yA(I=3) P(I=4) yW(I=5) y0(I=6) y1(I=7)
c     these are the data in Table 1 of Verner+ (1996,ApJ,465,487)
c     except eth is stored in PH1 for is=ntot(I) (see above)



*********************************************************************************

      implicit none
      integer           nz,ne,is,nout,nint
      double precision  s,e,einn,p1,x,z,y,q,a,b 


      integer           l,ninn,ntot
      double precision  ph1,ph2
      common/l/         l(7)           ! angular momentum indexed by shell number
      common/ninn/      ninn(30)       ! inner shell 
      common/ntot/      ntot(30)       ! number of shells for number of electrons
      common/ph1/       ph1(6,30,30,7) ! fitting parameters, inner shells
      common/ph2/       ph2(7,30,30)   ! fitting parameters, outer shells


c     set the cross section to the default null value

      phfit2 = 0.0d0

c     if the species is off the table then return the null cross section

      IF (nz.lt.1.OR.nz.gt.30) RETURN

c     if the number of electrons is out of range for the species then
c     retrun the null cross section

      IF (ne.lt.1.OR.ne.gt.nz) RETURN

c     determine the number of the outer valence electron shell for this
c     ion; basically nout is the shell number of the outermost shell
c     that is populated by an electron for this ion (assuming ground
c     state)

      nout = ntot(ne)
      IF (nz.eq.ne.AND.nz.gt.18) nout = 7
      IF (nz.eq.(ne+1).AND.(nz.eq.20.OR.nz.eq.21.OR.nz.eq.22.OR.nz.
     &   eq.25.OR.nz.eq.26)) nout = 7

c     this accounts for the fact that the 3d (is=6) shell of KI, CaI,
c     and CaII do not fill (though the 4s shell does (is=7); the matrix
c     PH1 reflects this properly

      IF (nz.eq.19.AND.ne.eq.19.AND.is.eq.6) RETURN
      IF (nz.eq.20.AND.ne.eq.19.AND.is.eq.6) RETURN
      IF (nz.eq.20.AND.ne.eq.20.AND.is.eq.6) RETURN

c     if we input an electron shell=is greater than the outer most
c     populated electron shell for this ion, nout, then return the null
c     cross section; the maximum value of nout=7 and therefore the
c     maximum value of is=7

      IF (is.gt.nout) RETURN

c     if the energy is below the ionization threshold energy for this
c     shell of this ion, then return the null cross section; matrix
c     element parameter ph1(1,nz,ne,is) is the threshold energy for
c     ionization for shell=is of this ion; if is=ntot(I) (where I
c     denotes nz=k) then this is the outer valence electron threshold
c     'eth' as tabulated in Table 1.  if is=1, then this is the most
c     inner electron and this is 'emax' as tabulated in Table 1; for
c     is>1 and is<ntot(I), this is the threshold for the shell=is

      IF (e.lt.ph1(1,nz,ne,is)) RETURN

c     determine the energy, einn, of the inner most shell, nint; there
c     is some logical directives going on here with regard to the value
c     assigned to einn (for the final calculation, see below) 

c     (1) setting einn=0 occurs if we do not have entires in the PH2
c     table for the value of nz, then e>einn (see below) and use of the
c     PH1 table is enforced; nz=15,17,19,21-25,27-30 NOT in PH2 table

c     (2) if ne<3, then we are working with H or He, so we then set
c     einn=10^30 so that e<einn is assured- then the condition is<=nint
c     uses PH1 and the condition is>nint uses PH2 (basically, PH2 use
c     means we are in the outermost shell, which is the ground state
c     ionization

c     (3) lastly, we set einn=ionization energy of the inner most shell,
c     then if the choice depends upon either is<nint OR e>einn for PH1;
c     if both are false then we are in the outer shell (which is the
c     ground state ionization) and we use PH2

      nint = ninn(ne)

      IF (nz.eq.15.OR.nz.eq.17.OR.nz.eq.19.or. 
     & (nz.gt.20.AND.nz.ne.26)) then
       einn = 0.0d0                           
      ELSE
       IF (ne.lt.3) then       ! ensure nint>=1  for PH1 index
        einn = 1.0d+30         ! arbitrary large number (e threshold)
       ELSE
        einn=ph1(1,nz,ne,nint) ! eX for nz,ne,nint
       END IF
      END IF

c     AND gate: nint<is<nout electron shell between min and max
c               e<eX(nint)   energy condition (F for nz=15,17,19,21-15,27-20)

      IF (is.lt.nout.AND.is.gt.nint.AND.e.lt.einn) return

c     finally! we are on the grid and in the allowable energy regime; if
c     the shell is less than the inner valence shell OR the energy is
c     greater than the threshold of the inner most shell, then compute
c     the partial xsec (first block); otherwise compute the full xsec
c     (lower block)

      IF (is.le.nint.OR.e.ge.einn) then
       p1 = -ph1(5,nz,ne,is)                                ! -P
       y  = e/ph1(2,nz,ne,is)                               ! x = e/e0 (y0=0)
       q  = -0.5*p1-l(is)-5.5                               ! modified power of y
       a  = ph1(3,nz,ne,is)*((y-1.0)**2+ph1(6,nz,ne,is)**2) ! s0*F(y) (yW term)
       b  = sqrt(y/ph1(4,nz,ne,is))+1.0                     ! F(y)    (yA term)
       s  = a*y**q*b**p1                                    ! full s0*F(y)
      ELSE
       p1 = -ph2(4,nz,ne)                                   ! -P
       q  = -0.5*p1-5.5                                     ! power of y
       x  = e/ph2(1,nz,ne)-ph2(6,nz,ne)                     ! x = e/e0 - y0
       z  = sqrt(x*x+ph2(7,nz,ne)**2)                       ! y = sqrt(x^2+y1^2)
       a  = ph2(2,nz,ne)*((x-1.0)**2+ph2(5,nz,ne)**2)       ! s0*F(y) (yW term)
       b  = 1.0+sqrt(z/ph2(3,nz,ne))                        ! F(y)    (yA term)
       s  = a*z**q*b**p1                                    ! full s0*F(y)
      END IF

c     set the cross section = to the function 

      phfit2 = s

c     return

      RETURN

      END
