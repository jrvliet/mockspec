c     rates-alphas.f
c     routines to populated the rate coefficients

c     included routines:
c     SUBROUTINE alphas_nonT(k,j)
c     SUBROUTINE alphas_T(k,j,T)
c
c.........................................................................
c

      SUBROUTINE alphas_nonT

c     obtain the alphas (rate coefficients) for ion k,j that are
c     not temperature dependent

c     this routine called only if doPH=.true.
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           j,k,kk
      double precision  PH_Rate,AUG_Rate
      include           'rates.com'
      include           'com-modes.com'



c     loop over the species and ionization stages and compute the rate
c     coefficients for all included processes; loop over j ionization
c     stages for each species k, we do not make rates for the fully
c     ionized ion; so we go to the kth ionization levels (note that the
c     ionization rates are indexed by the initial ion)


      DO kk=1,Nspecies
      
       k = kidx(kk)

       DO j=1,k

c     ...............
c     PHOTOIONIZATION
c     ...............

c     compute the photoionization rate coefficients (function PH_Rate);
c     the photoionization rate requires integrating over the ionizing
c     spectrum and thus is dependent upon the ionizing photon density;
c     the photoionization rate is weighted by the Auger probability
c     yield for single electron ejection; function PH_Rate also stores
c     the partial unweighted photoionization rates for each shell of ion
c     k,j; these unweighted partial rates are then used to compute the
c     Auger rates

        R_ph(k,j) = PH_Rate(k,j)

c     ................
c     AUGER IONIZATION
c     ................

c     compute the Auger rates (subroutine getQs); as with the
c     photoionization rates, the Auger rates require integrating over
c     the ionizing spectrum and are thus also dependent upon the
c     ionizing photon density; the partial photoionization rates
c     returned by function PH_Rate are weighted by the Auger probability
c     yields for various multiple electron ejections; finally, for
c     simplification of Auger destruction rates, we compute and the
c     summed rate coefficient for all Auger ioniations out of k,j
c     (function AUG_rate)

        IF (doAUG) then
          CALL getQs(k,j)
          R_Agrout(k,j) = AUG_Rate(k,j)
        ENDIF

       ENDDO  ! next j

      ENDDO   ! next k


      RETURN

      END


c
c.........................................................................
c

      SUBROUTINE alphas_T(T)

c     obtain the alphas (rate coefficients) for ion k,j that are
c     temperature dependent

c     IMPORTANT CONVENTION

c     all ionization rate coefficients are indexed by k,j for
c     ionization from k,j to k,j+1 (by the initial ion)

c     all recombination rate coefficients are indexed by k,j for
c     recombination from k,j+1 to k,j (by the final ion)

c     thus there are no rates indexed for the fully ionized ion of
c     species k, in which j=k+1
c     

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           j,k,kk
      double precision  T
      double precision  REC_Rate,CDI_Rate,CEA_Rate
      double precision  HCTRecomb,HCTIon,HECTRecomb
      include           'rates.com'
      include           'com-modes.com'


c     loop over the species and ionization stages and compute the rate
c     coefficients for all included processes; loop over j ionization
c     stages for each species k, we do not make rates for the fully
c     ionized ion; so we go to the kth ionization levels (note that the
c     recombination rates are indexed by the final ion, not the initial
c     ion, whereas the ionization rates are indexed by the initial ion)


      DO kk=1,Nspecies
      
       k = kidx(kk)

       DO j=1,k

c     .............
c     RECOMBINATION
c     .............

c     compute the TOTAL recombination rate coefficient (function
c     REC_Rate); this rate includes photo, low T dielectronic, and high
c     T dielectronic recombinations

        alpha_rec(k,j) = REC_Rate(k,j,T)

c     ......................
c     COLLISIONAL IONIZATION
c     ......................

c     compute the direct collision ionization rate (function CDI_Rate)
c     and the excitation-autoionization collisional ionization rate
c     (function CEA_Rate)

        IF (doCDI) then
         alpha_Cdi(k,j) = CDI_Rate(k,j,T)
         IF (doCEA) alpha_Cea(k,j) = CEA_Rate(k,j,T)
        ENDIF

c     ...............
c     CHARGE TRANSFER
c     ...............

c     compute the charge transfer collisional ionization and
c     recombination rates

        IF (doCT) then
         alpha_CTrecH(k,j)  = HCTRecomb(k,j,T)
         alpha_CTionH(k,j)  = HCTIon(k,j,T)
         alpha_CTrecHe(k,j) = HeCTRecomb(k,j,T)
        ENDIF

       ENDDO  ! next j

      ENDDO   ! next k



      RETURN
      END
