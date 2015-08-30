c     PHOTOIONIZATION ARRAYS FOR IONIZING SPECTRA AND SELF-SHIELDING

c     SIZES (LOGICAL) OF ARRAYS
c     NExsec = number of array elements for the cross section energy
c     NEuvb  = number of array elements for the UVB spectrum energy
c     NEsb99 = number of array elements for the Starburst99 spectrum energy
c     NEsed  = number of array elements for the employed spectrum energy

      integer           NExsec,NEuvb,NEsb99,NEsed,kp,jp,sp
      COMMON/intblk/    NExsec,NEuvb,NEsb99,NEsed,kp,jp,sp


c     RANGE OF PHOTON ENERGIES (for integrations to obtain rates)
c     logEmin = minimum energy (defined by ion k,j threshold energy)
c     logEmax = maximum energy (defined by SED maximum energy)
c     Ebuff   = buffer for interpolation purposes

      double precision  logEmin,logEmax,Ebuff
      COMMON/Eblk/      logEmin,logEmax,Ebuff


c     UVB DATA (Haardt & Madau 1996 ultraviolet background)
c     zuvb      = array of grid z for logJnu array matrix element
c     logEuvb   = array of grid log(Euvb) for logJNEuvb matrix element
c     logJNUuvb = log(Jnu) matrix of redshift z and energy E grid
c     logJEuvb  = log(JE) array at grid log(Euvb) interpolated for desired z


      double precision  zuvb,logEuvb,logJNUuvb,logJEuvb
      COMMON/UVBblk/    zuvb(UVBzpix),logEuvb(UVBmaxpix),
     &                  logJNUuvb(UVBzpix,UVBmaxpix),
     &                  logJEuvb(UVBmaxpix)


c     SB99 DATA (Starburst 99 stellar population spectra)
c     Zsolsb99    = array of metallicities for grid stellar population
c     Agesb99     = array of ages for SB99 grid stellar population
c     logEsb99    = array of grid log(Esb99) for SB SED element
c     logJEsb99   = log(JE) array at grid log(Euvb) for target stellar pop
c     logLlamsb99 = array of log(Llamsb99), input luminosity density array 

      double precision  Zsolsb99,Agesb99,Msolsb99,logEsb99,
     &                  logJEsb99,logLlamsb99
      COMMON/SB99blk/   Zsolsb99(SB99NZ),
     &                  Agesb99(SB99Nage),
     &                  logEsb99(SB99maxpix),
     &                  logJEsb99(SB99maxpix),
     &                  logLlamsb99(SB99Nage,SB99NZ,SB99maxpix)


c     SED DATA (adopted SED used for the photoionization calcs)
c     logJEsed  = log(Jnu) array energy E grid (for all cals)
c     logEsed   = array of grid log(Esed) for logJEsed array element
c     d2JEsed   = 2nd derivitives for spline interpolation

      double precision  logJEsed,logEsed,d2JEsed
      COMMON/SEDblk/    logJEsed(SEDmaxpix),
     &                  logEsed(SEDmaxpix),
     &                  d2JEsed(SEDmaxpix)


c     VARIABLES FOR OPTICALLY THICK TREATMENT
c     E11,sig11 = neutral hydrogen threshold energy and cross section
c     E21,sig21 = neutral helium threshold energy and cross section
c     E22,sig22 = singly ionized helium threshold energy and cross section

      double precision  E11,sig11,E21,sig21,E22,sig22,taumin
      COMMON/xsecblk/   E11,sig11,E21,sig21,E22,sig22,taumin
