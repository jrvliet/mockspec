c
c.........................................................................
c

      SUBROUTINE getQs(k,j)

c
c     author: Chris Churchill
c
c     compute the Auger Q coefficients
c
c     see Eq 10 from Verner + Yakovlev (1990,A&SS,165,27)
c
c
c     k   = the species number
c     j   = the ionization stage of the pre ionized ion
c     m   = the ionization stage of the post ionized ion
c     s   = the shell number
c     ne  = the number of electrons ejected going from j to m
c
c     m is always greater than j+2 because j+1 yields ne=1, and ne=1 is
c     the non-Auger single electron photoionization, which is stored
c     globally as R_ph and was computed in function PH_Rate
c
c     the partial photoionization rate for shell "s", R_phs, was
c     also computed in function PH_Rate
c
c     the quantity W_Auger is the probability that ne electrons are
c     ejected from shell "s" from initial ion k,j (which subsequently
c     becomes ion k,m)
       
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           k,j,m,s,ne,nek,nout
      double precision  asum
      include           'rates.com'
      include           'com-auger.com'

c     we need the total number of populated shells, so include the ntot
c     common block

      integer           ntot
      COMMON/ntot/      ntot(30)



c     for the initial ion k,j, determine the number of bound electrons,
c     nek, and the index of the most outer shell populated by one or
c     more electron(s) for this number of electrons, nout; nout is
c     determined from Verner's ntot array; it is assumed the ion is in
c     the ground state

      nek  = k - j + 1
      nout = ntot(nek)

c     loop over the final ionization stages (m>=j+2) following the Auger
c     process from initial ionization stage j; compute the number of
c     ejected electrons; null the partial rate running sum (asum), loop
c     over the populated shells of the initial ion (j) and sum the
c     partial rates to obtain the total rate coefficient, Q

      DO 11 m=j+2,k+1
       ne    = m - j
       asum  = 0.0d0
       IF (ne.le.10) then  ! tables stop at 10 ejected electrons
        DO 15 s=1,nout
         asum = asum + W_Auger(k,j,s,ne)*R_phs(s)
 15     CONTINUE
       END IF
       Q(k,j,m) = asum
 11   CONTINUE


      END


c
c.........................................................................
c

      DOUBLE PRECISION FUNCTION AUG_Rate(k,j)

c
c     author: Chris Churchill
c
c
c     compute the Auger rate for destrcution of k,j
c
c     k   = the species number
c     j   = the ionization stage of the pre ionized ion
c     m   = the ionization stage of the post ionized ion
c
c     m is always greater than j+2 because j+1 yields ne=1, and ne=1 is
c     the non-Auger single electron photoionization, which is stored
c     globally as R_ph and was computed in function PH_Rate

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           k,j,m
      double precision  asum
      include           'rates.com'


      asum =  0.0d0

c     sum over the final ionization stages m out of which initial
c     ionization stage j can go to following Auger ionization of 2 or
c     more electrons (where the photoionized electron is the first
c     electron)

      DO 11 m=j+2,k+1
       asum = asum + Q(k,j,m)
 11   CONTINUE

      AUG_Rate = asum

      RETURN
      END
