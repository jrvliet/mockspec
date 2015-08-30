c     AUGER IONIZATION YIELD PROBABILITIES

c     the yeild probabilities are from Kaastra + Mewe (1993.A&ASS,97,443)
c     for ionization of Nej ejected electrons from shell s of ion k,j
c     the tabulated data range from Nej=1 to Nej=10 

c     W_Auger = yield probabilities for k,j,s,Nej

      double precision  W_Auger
      COMMON /rauger/   W_Auger(Imax,Imax,Nshmax,10)
     &                  

