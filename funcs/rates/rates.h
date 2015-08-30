c     PARAMTER FILE FOR RATES.F



c     SPECIES/ION GRID
c............................................................................
c     Imax  = maximum number of ions
c     NSHmax = maximum number of shells (Verners)
c     STGmax = maximum number of ionization stages for recombination 
c     Smax   = maximum number of shells for direct collisional ionization 

      integer    Imax,Nmax,NSHmax,STGmax,Smax
      parameter  (Imax   =     30,
     &            Nmax   =    496,
     &            NShmax =      7,
     &            STGmax = Imax-1,
     &            Smax   =      3)



c     UVB GRID
c............................................................................
c     UVBmaxpix = number of lines in Haardt & Madau 1996 SEDs
c     UVBzpix  = number of H&M96 files (each one at a given z)

      integer    UVBmaxpix,UVBzpix
      parameter  (UVBmaxpix = 3750,
     &            UVBzpix   =   26)



c     SB99 GRID
c............................................................................
c     SB99maxpix = number of lines in Starburst99 spectra SEDs
c     SB99Nage   = number of ages of SB99 SEDS
c     SB99NZ     = number of metallicities of SB99 SEDS
c     SB99maxAZ  = maximum of SB99Nage and SB99NZ (tmp storage array)

      integer    SB99maxpix,SB99Nage,SB99NZ,SB99maxAZ
      parameter  (SB99maxpix =     1221,
     &            SB99Nage   =        5,
     &            SB99NZ     =        2,
     &            SB99maxAZ  = SB99Nage)



c     SED IONIZING SPECTRUM GRID
c............................................................................
c     SB99maxpix = number of lines in Starburst99 spectra SEDs
c     SB99types  = number of SB99 SEDS

      integer    SEDmaxpix
      parameter  (SEDmaxpix = UVBmaxpix)



c     ENERGY AND TEMPERATURE GRIDS
c............................................................................
c     NEmax - maximum number of energy points on E grid

      integer    NEmax
      parameter  (NEmax = 9000)


c     PHYSICAL CONSTANTS AND CONVERSIONS
c............................................................................
c     h      = Planck constant  (erg sec)
c     clight = speed of light (cm/s)
c     pi     = pie, yummm
c     erg2eV = conversion of 1 eV to ergs
c     pc2cm  = conversion of 1 pc to cm

      double precision  h,clight,pi,erg2eV,pc2cm
      parameter         (h      =       6.626075540e-27,       
     &                   clight =         2.99792458d10, 
     &                   pi     = 3.1415926535897932384, 
     &                   erg2eV =        1.60217646d-12,
     &                   pc2cm  =         3.08568025d18) 
