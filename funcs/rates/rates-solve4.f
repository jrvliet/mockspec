c     rates-solve4
c     routines to solve the ionization balance of the gas

c     included routines:
c     SUBROUTINE ionbalance(toler)
c     SUBROUTINE ionbal0(ne0,toler,error)
c     FUNCTION   zerometals(ne)
c     FUNCTION   econserve(ne)
c     SUBROUTINE solve_metals(ne)
c     SUBROUTINE refine_fHfHe(ne)

c
c.........................................................................
c

      SUBROUTINE ionbalance(toler)

c     this routine is the solver for the ionization balance by zeroing a
c     series of linearized matrix equation using the principle of charge
c     conservation as the constrain; the zeroing (charge conservation)
c     is performed to an accuracy=toler, which is passed to this routine

c     STEP 1. obtain a zeroth order ionization solution for the electron
c     density and ionization fraction of the hydrogen and helium ions
c     (omiting metals); using these values, obtain a zeroth order
c     solution of the ionization fractions of all metal ions; this is
c     the "initial guess" model on which we iterate toward solution

c     STEP 2. using charge conservation of free electrons, solve for the
c     equilibrium ionization fractions for all ions; then compute the
c     equilibrium number densities of all metal ions; return

c     Note: for ions k,j that have no recombination rates, the program
c     will crash!

c     the ionization fractions are global variables

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      external          econserve
      logical           error
      double precision  econserve,zbrent
      integer           kk,k,j
      double precision  toler
      double precision  ne0,ne_min,ne_max
      include           'rates.com'


c     initialize error flag, which is set high when the solver fails to
c     converge

      error = .false.

c     ........................
c     STEP 1: zeroth solution
c     ........................

c     obtain the zeroth estimate for the electron density and ionization
c     fractions of the hydrogen and helium ions; on return, if
c     error=.true., notify and abort

      CALL ionbal0(ne0,toler,error)

      IF (error) then
        WRITE(6,*) 'ERROR generated in subroutine ionbal0'
        WRITE(6,*) 'while zeroing function zerometals'
        WRITE(6,*) 'program aborted'
        WRITE(4,*) 'ERROR generated in subroutin ionbal0'
        WRITE(4,*) 'while zeroing function zerometals'
        WRITE(4,*) 'program aborted'
        STOP 
      ENDIF

c     solve for the zeroth order ionization fractions of the metal ions
c     using charge conservation assuming the 0th solution for the
c     electron density and hydrogen and helium ionization fractions

      CALL solve_metals(ne0)

c     .............................
c     STEP 2: equilibrium solution
c     .............................

c     we are now armed with the initial starting solution

c     first, define the range over which to search for the equilibrium
c     value of the electron density; the minumum is conservatively set,
c     the maximum is for a fully ionized gas with the input metallicity

      ne_min = 0.10d0*ne0     

      ne_max = 0.0d0
      DO kk=1,Nspecies
       k = kidx(kk)
       ne_max = ne_max + real(k)*nk(k) 
      ENDDO
      ne_max = 10.0d0*ne_max

c     adopting the zeroth solution as the starting location, use Brent's
c     method to iterate toward the equilibrium solution by adjusting
c     until charge density is conserved to an accuracy=toler; the
c     function "econserve" evaulates the charge conservation; on return,
c     if error=.true., notify and abort

      eden = zbrent(econserve,ne_min,ne_max,toler,error)

      IF (error) then
        WRITE(6,*) 'ERROR generated in subroutine ionbalance'
        WRITE(6,*) 'while zeroing function econserve'
        WRITE(6,*) 'program aborted'
        WRITE(4,*) 'ERROR generated in subroutin ionbalance'
        WRITE(4,*) 'while zeroing function econserve'
        WRITE(4,*) 'program aborted'
        STOP 
      ENDIF

c     from the equilibrium ionization fractions, compute the number
c     densities of all ions 

      DO kk=1,Nspecies
       k = kidx(kk)
       DO j=1,k+1
         nkj(k,j) = fion(k,j)*nk(k)
       ENDDO
      ENDDO

c     that was easy!

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE ionbal0(ne0,toler,error)

c     called from SUBROUTINE ionbalance
c     employs FUNCTION zbrent, which evaluates FUNCTION zerometals

c     this routine is the driver for obtaining the zeroth order solution
c     to the ionization balance assuming a gas comprising hydrogen and
c     helium only; the solution is found employing charge conservation

c     we return "ne0", the zeroth order estimate of the electron density
c     and ionization fractions the hydrogen and helium

c     the ionization fractions are global variables

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      external          zerometals
      logical           error
      double precision  zerometal,zbrent
      double precision  ne0,ne_min,ne_max,toler
      include           'rates.com'


c     bracket the 0th order electron density; the minimum, ne_min, is ne
c     for a virtually neutral gas and the maximum,ne_max, is ne for a
c     fully ionized (again, no metals)

      ne_min = 1.0d-15*nk(1)             ! ne for virtually neutral
      ne_max = 1.0d0*nk(1) + 2.0d0*nk(2) ! ne for fully ionized 

c     apply charge conservation to obtain the ne0 and the zeroth order
c     ionization fractions for hydrogen and helium

      ne0 = zbrent(zerometals,ne_min,ne_max,toler,error)

      RETURN
      END

c
c.........................................................................
c

      DOUBLE PRECISION FUNCTION zerometals(ne)

c     called by FUNCTION zbrent in SUBROUTINE ionbal0
c     returns the zeroth order electron density 

c     this function is to be zeroed (root solved) for the zeroth order
c     electron density ne=ne0 for a zero metallicity gas (hydrogen and
c     helium only); we also obtain the zeroth order estimate of
c     ionization fractions for hydrogen and helium; this solution does
c     NOT include charge exchange with metals
c
c     for the adjacent number density ratios we adopt the notation:
c
c     for hydrogen: phik1=n(1,2)/n(1,1)
c     for helium:   phik1=n(2,2)/n(2,1) , phik2=n(2,3)/n(2,2)

c     the ionization fractions are global variables

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           j,k
      double precision  phik1,phik2,ne,nesum
      include           'rates.com'

c     hydrogen: adjacent ionization stage ratio of number densities and
c     ionization fractions; no charge exchange

      phik1   = (R_ph(1,1) + ne*(alpha_Cdi(1,1)+alpha_Cea(1,1)))
     &         /(ne*alpha_rec(1,1))

      fion(1,1) = 1.0d0/(1.0d0+phik1)
      fion(1,2) = phik1*fion(1,1)

c     helium: adjacent ionization stage ratio of number densities and
c     ionization fractions; no charge exchange

      phik1   = (R_ph(2,1) + ne*(alpha_Cdi(2,1)+alpha_Cea(2,1)))
     &         /(ne*alpha_rec(2,1))
      phik2   = (R_ph(2,2) + ne*(alpha_Cdi(2,2)+alpha_Cea(2,2)))
     &         /(ne*alpha_rec(2,2))

      fion(2,1) = 1.0d0/(1.0d0+phik1+phik1*phik2)
      fion(2,2) = phik1*fion(2,1)
      fion(2,3) = phik2*fion(2,2)
      
c     determine the electron density using charge conservation

      nesum = 0.0d0
      DO 11 k=1,2     ! sum only over hydrogen and helium
       DO 13 j=2,k+1  ! neutral does not contribute electrons
        nesum = nesum + real(j-1)*fion(k,j)*nk(k)
 13    CONTINUE
 11   CONTINUE

c     does ne=nesum? (root solve); perform in log space for more stable
c     convergence by routine zbrent

      zerometals = log10(ne) - log10(nesum) 

      RETURN
      END


c
c.........................................................................
c

      DOUBLE PRECISION FUNCTION econserve(ne)

c     evaluated by FUNCTION zbrent, which is called by SUBROUTINE ionbalance
c     calls SUBROUTINE refine_fHfHe and SUBROUTINE solve_metals
c     each iteration returns the current estimate of the electron density

c     THE MAIN SOLVER: conserve charge density to solve the ionization
c     fractions of all ions

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           k,kk,j
      double precision  ne,nesum
      include           'rates.com'


c     refine hydrogen and helium ionization fractions by accounting for
c     charge exchange with metals

      CALL refine_fHfHe(ne)

c     solve for ionization fractions of the metals, including charge
c     exchange with hydrogen and helium and Auger ionizations

      CALL solve_metals(ne)

c     obtain the summed electron density by employing charge
c     conservation

      nesum = 0.0d0
      DO kk=1,Nspecies
       k = kidx(kk)
       DO j=2,k+1  ! the neutral stage does not contribute
        nesum = nesum + real(j-1)*fion(k,j)*nk(k)
       ENDDO
      ENDDO

c     is charge density conserved?  use logs for increased stability

      econserve = log10(ne) - log10(nesum)

      RETURN
      END


c........................................................................
c
      SUBROUTINE solve_metals(ne)

c     called by FUNCTION econserve
c     returns ionization fractions of all metals

c     this routine solves for the current estimate of the metal
c     ionization fractions using the current estimates of the electron
c     density and the hydrogen and helium ionization fractions; the
c     routine is broken into two steps

c     STEP 1: populate the rate matrix for all metal ions

c     STEP 2: using linear backsubstitution, we compute the phik(j) for
c     all ions, where phik(j) is the ratio of the number density of
c     adjacent ionization stages (j+1 over j) for species k; as we
c     proceed through each species k, we compute the ionization
c     fractions of all ionization stages j

c     the ionization fractions are global variables
c     the phik(j) are local variables

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit none

      include           'rates.h'
      integer           k,kk,j,i
      double precision  ne
      double precision  R_auger,x,y
      double precision  phik(Imax),Sk,Pk(Imax)
      include           'rates.com'


c     rate matrix notation: (global variables)
c     R_ion(k,j) = destruction of k,j to k,j+1 via ionization
c     R_rec(k,j) = creation of k,j via recombination to k,j-1
c     R_hat(k,j) = total destruction of k,j by all processes

c     .................................
c     STEP 1: populated the rate matrix
c     .................................

      DO kk=3,Nspecies
       k = kidx(kk)

c     j=1 (neutral ion) is a special case

       j = 1
       R_ion(k,j) = R_ph(k,j)                            
     &             + ne*(alpha_Cdi(k,j)+alpha_Cea(k,j))  
     &             + fion(1,2)*nk(1)*alpha_CTionH(k,j)   
       R_rec(k,j) = 0.0d0  ! never accessed              
       R_hat(k,j) = R_ion(k,j) + R_Agrout(k,j)           

c     all cases from j=2 to j=k

       DO j=2,k
        R_ion(k,j) = R_ph(k,j) 
     &              + ne*(alpha_Cdi(k,j)+alpha_Cea(k,j)) 
     &              + fion(1,2)*nk(1)*alpha_CTionH(k,j)
        R_rec(k,j) = ne*alpha_rec(k,j-1) 
     &              + fion(1,1)*nk(1)*alpha_CTrecH(k,j-1) 
     &              + fion(2,1)*nk(2)*alpha_CTrecHe(k,j-1)
        R_hat(k,j) = R_ion(k,j) + R_rec(k,j) + R_Agrout(k,j)
       ENDDO

c     j=k+1 (fully ionized ion) is a special case

       j = k + 1
       R_ion(k,j) = 0.0d0  ! never accessed
       R_rec(k,j) = ne*alpha_rec(k,j-1) 
     &             + fion(1,1)*nk(1)*alpha_CTrecH(k,j-1) 
     &             + fion(2,1)*nk(2)*alpha_CTrecHe(k,j-1)
       R_hat(k,j) = 0.0d0  ! never accessedC


      ENDDO  ! next k


c     ......................................................
c     STEP 2: solve for the phik(j) and ionization fractions
c     ......................................................


c     the loop labeled "Auger ionization contributions" computes the
c     summed rate of all Auger ionization of k,j (destruction term)
c     this sum is stored as R_auger

      DO kk=3,Nspecies

       k = kidx(kk)

c     j=1 and j=2: trivial cases, solve manually

       phik(1) = R_hat(k,1)/R_rec(k,2)
       phik(2) = (R_hat(k,2)-(R_ion(k,1)/phik(1)))/R_rec(k,3)

c     j=3 to j=k [there is no phik(j=k+1)]

       DO j=3,k       
        R_auger = 0.0d0
        x = phik(j-1)
        DO i=j-2,1,-1   ! Auger ionization contributions
         x   = phik(i)*x
         R_auger = R_auger + Q(k,i,j)/abs(x)
        ENDDO
        y       = (R_ion(k,j-1)/phik(j-1))/R_hat(k,j)+R_auger/R_hat(k,j)
        phik(j) = (R_hat(k,j)/R_rec(k,j+1))*y*(1.0d0/y-1.0d0)
C        if (phik(j).lt.0.) then
C         write (99,*) ionID(k,j)(1:6),phik(j)
C         phik(j) = 1.0d-30
C        endif    
       ENDDO

c     compute the denominator (Sk) of the ionization fractions

       Sk    = 0.0d0
       Pk(1) = 1.0d0
       DO j=2,k+1
        Pk(j) = Pk(j-1)*abs(phik(j-1))
       ENDDO      
       DO j=1,k+1
        Sk = Sk + Pk(j)
       ENDDO

c     compute the ionization fractions

       fion(k,1) = 1.0d0/Sk
       DO j=2,k+1
        fion(k,j) = abs(phik(j-1))*fion(k,j-1)
       ENDDO


      ENDDO  ! next k
       

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE refine_fHfHe(ne)

c     called by FUNCTION econserve
c     returns ionization fractions of hydrogen and helium

c     this routine solves for the hydrogen and helium ionization
c     fractions for the current iteration toward the solution using the
c     current estimate of the electron density and the previous
c     iteration's metal ionization fractions; that is, the refinement is
c     the updating of the charge exchange with metals based upon the
c     previous iteration of the equilibrium solution

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'      
      integer           kk,k,j
      double precision  ne
      double precision  top,bot,phi1,phi2
      double precision  RexH_11,RexHp_12,RexHe_21,RexHp_22,RexH_23
      include           'rates.com'


c     notation for hydrogen balance charge exchange rates
c     RexH_11  = ionization destruction of neutral hydrogen by charge
c                exchange recombination to metal ion k,j (R^{exH}_{1,1})
c     RexHp_12 = recombination creation of neutral hydrogen by charge
c                exchange ionization of metal ion k,j (R^{exH+}_{1,2})

c     notation for helium balance charge exchange rates
c     RexHe_21 = ionization destruction of neutral helium by charge 
c                exchange recombination to metal ion k,j (R^{exHe}_{2,1})
c     RexHp_22 = creation of singly ionized helium by charge exchange
c                recombination with ionized hydrogen (R^{exH+}_{2,2})
c     RexH_23  = destruction of doubly ionized helium by charge exchange
c                ionization of neutral hydrogen (R^{exH}_{2,3})

      RexH_11  = 0.0d0
      RexHp_12 = 0.0d0
      RexHe_21 = 0.0d0  
      RexHp_22 = 0.0d0
      RexH_23  = 0.0d0

c     obtain the hydrogen balance charge exchange rates summed over
c     helium and all metals

      DO kk=2,Nspecies
       k = kidx(kk)
       DO j=2,k+1
        RexH_11 = RexH_11 + fion(k,j)*nk(k)*alpha_CTrecH(k,j-1)  
       ENDDO
       DO j=1,k
        RexHp_12 = RexHp_12 + fion(k,j)*nk(k)*alpha_CTionH(k,j)
       ENDDO
      ENDDO

c     obtain the helium balance charge exchange rate summed over all
c     metals; we must first obtain the hydrogen balance before we can
c     add/compute charge exchange with hydrogen (see below)

      DO kk=3,Nspecies
       k = kidx(kk)
       DO j=2,k+1
        RexHe_21 = RexHe_21 + fion(k,j)*nk(k)*alpha_CTrecHe(k,j-1)
       ENDDO
      ENDDO

c     solve hydrogen

      top  = R_ph(1,1) + ne*alpha_Cdi(1,1) + RexH_11
      bot  = ne*alpha_rec(1,1) + RexHp_12
      phi1 = top/bot

      fion(1,1) = 1.0d0/(1.0d0+phi1)
      fion(1,2) = phi1*fion(1,1)

c     a sanity check

      IF (fion(1,1).gt.1.0) then
       WRITE(6,*) 'ERROR(refine_fHfHe): unphysical fion(1,1)'
       WRITE(6,*) 'fion11 = ',fion(1,1)
       WRITE(6,*) 'phi11  = ',phi1
       STOP
      END IF

c     solve helium; we first need to add hydrogen charge exchange to
c     RexHe_21 and compute RexHp_22 and RexH_23

      RexHe_21 = RexHe_21 + fion(1,2)*nk(1)*alpha_CTrecHe(1,1) 
      RexHp_22 = fion(1,2)*nk(1)*alpha_CTionH(2,2)
      RexH_23  = fion(1,1)*nk(1)*alpha_CTrecH(2,2)

      top  = R_ph(2,1) + ne*alpha_Cdi(2,1) + RexHe_21
      bot  = ne*alpha_rec(2,1)
      phi1 = top/bot

      top  = R_ph(2,2) + ne*alpha_Cdi(2,2) + RexHp_22
      bot  = ne*alpha_rec(2,2) + RexH_23
      phi2 = top/bot

      fion(2,1) = 1.0d0/(1.0d0+phi1+phi1*phi2)
      fion(2,2) = phi1*fion(2,1)
      fion(2,3) = phi2*fion(2,2)


      RETURN

      END

c     eof
c
c.........................................................................
c
