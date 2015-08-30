c     GRIDS AND ARRAYS FOR SELF-SHIELDING 

c     GRID FOR SOLVING SELF SHIELDING
c     for each k,j,s there is a 3D grid as a function of hydrogen density,
c     temperature, and incremental optical depth at the dominant ionization
c     edge from HI, HeI, or HeII

c     NRnmax   = physical size of hydrogen number density, nH, axis
c     NRTmax   = physical size of temperature, T, axis
c     NRtaumax = physical size of optical depth (tau) axis

      integer           NRnmax,NRTmax,NRtaumax
      parameter         (NRnmax   = 50,
     &                   NRTmax   = 50, 
     &                   NRtaumax = 30)
      
c     NRn    = number of array elements for nH grid axis
c     NRT    = number of array elements for T grid axis
c     NRtau  = number of array elements for tau grid direction
c     nR     = array of hydrogen densities, nH
c     TR     = array of temperatures, T
c     tau    = array of incremental optical depths

      integer           NRn,NRT,NRtau
      COMMON/iRgrid/    NRn,NRT,NRtau      

      double precision  nR,TR,tau
      COMMON/RnTtaublk/ nR(NRnmax),TR(NRTmax),tau(NRtaumax)


c     THE GRID OF WEIGHTED MEAN IONIZATION FRACTIONS
c     fgrid = grid of ionization fractions for m,nH,T, where
c             m is the index of the target ions

      double precision  fgrid
      COMMON/fionblk/  fgrid(Imax,Nrnmax,NRTmax)


c     UNATTENUATED RATES (for optically thin cells, tau<1)
c     Rphkjs0 = 3D array of photoionization rates for shell s of ion k,j 
c               for tau<1

      double precision  Rphkjs0
      COMMON/Rzeroblk/  Rphkjs0(Imax,Imax,NShmax)
