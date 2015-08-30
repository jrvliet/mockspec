c     COMMON/GLOBAL VARIABLES FOR RATES.F


c     SPECIES, ABUNDANCES, IONIZATION FRACTIONS, ETC.

c     specflag = logical flag for incluse/exclusions of species k

      logical           specflag
      COMMON/logiblk/   specflag(Imax)

c     Nspecies = the number of species being model
c     Nprint   = number of cells processed between screen updates
c     kidx     = stores k index number for included atomic species
c     noutfiles = the number of ions to be output to grids

      integer           Nspecies,Nprint,noutfiles,kidx
      COMMON/ispecblk/  Nspecies,Nprint,noutfiles,kidx(Imax)

c     nH       = hydrogen number density
c     alphksol = abundance fraction of species k in solar mixture
c     fabund   = abundance fraction of species k in included mixture
c     xfrac    = mass fractions for species k in included mixture
c     nk       = number density of species k (all ionization stages)
c     nkj      = number density of ionization stage j of species k 
c     eden     = electron density
c     nions    = number density of all ions
c     fion     = ionization fraction of ion k,j

      double precision  nH,alphksol,fabund,xfrac,nk,nkj,eden,nions,fion
      COMMON/ions/      nH,alphksol(Imax),fabund(Imax),xfrac(Imax),
     &                  nk(Imax),nkj(Imax,Imax+1),eden,nions,
     &                  fion(Imax,Imax+1)


c     INDIVIDUAL OR PARTIAL RATE COEFFICIENTS 

c     R_phs         = photoionization rate for shell s
c     R_ph          = photoionization rate (ejecting 1 electron)
c     alpha_rec     = recombination rate coefficient (rad,loTde,hiTde)
c     alpha_Cdi     = direct collisional ionization rate coefficient
c     alpha_Cea     = ex-autoionization rate coefficient
c     alpha_CTrecH  = charge transfer H recombination rate coefficient
c     alpha_CTionH  = charge transfer H ionization rate coefficient
c     alpha_CTrecHe = charge transfer He recombination rate coefficent
c     R_Agrout      = total Auger ionization rate out of k,j
c     Q             = Auger ionization k,j ejecting >1 electron 

      double precision  R_phs        
      double precision  R_ph        
      double precision  alpha_rec,alpha_rec_ph,alpha_rec_die       
      double precision  alpha_Cdi       
      double precision  alpha_Cea       
      double precision  alpha_CTrecH    
      double precision  alpha_CTionH    
      double precision  alpha_CTrecHe   
      double precision  R_Agrout 
      double precision  Q      
      COMMON/ratblock/  R_phs(NShmax),
     &                  R_ph(Imax,Imax+1),
     &                  alpha_rec(Imax,Imax+1),
     &                  alpha_rec_ph(Imax,Imax+1),
     &                  alpha_rec_die(Imax,Imax+1),
     &                  alpha_Cdi(Imax,Imax+1),
     &                  alpha_Cea(Imax,Imax+1),
     &                  alpha_CTrecH(Imax,Imax+1),
     &                  alpha_CTionH(Imax,Imax+1),
     &                  alpha_CTrecHe(Imax,Imax+1),
     &                  R_Agrout(Imax,Imax+1),
     &                  Q(Imax,Imax+1,Imax+1)


c     TIME SCALES (stored by target ion index, not k,j ; see rates-times.f)

c     tau_ph   = photoionization timescale
c     tau_rec  = recombination timescale
c     tau_coll = collisional iomnization timescale
c     tau_cool = coronal approx (collisional only) cooling time

      double precision  tau_ph,tau_rec,tau_coll,tau_cool
      COMMON/taublk/    tau_ph(Imax),tau_rec(Imax),tau_coll(Imax),
     &                  tau_cool



c     RATES FOR SOLVER (see rates-solve4.f)

c     R_hat(k,j) = total destruction rate of k,j
c     R_ion(k,j) = total creation rate of k,j-1 by ionization of k,j
c     R_rec(k,j) = total creation rate of k,j+1 by recombination to k,j

      double precision  R_hat,R_ion,R_rec
      COMMON/totrats/   R_hat(Imax,Imax+1),
     &                  R_ion(Imax,Imax+1),
     &                  R_rec(Imax,Imax+1)


c     CHARACTER DATA FOR COMMUNICATIONS 

c     specID   = name of elemental species k,(i.e., k=12 is Mg)
c     ionID    = name of ion k,j (i.e., k=12,j=2 is MgII)
c     shellID  = shell spectroscopic notation (2s, etc)
c     group    = Group type on Periodic Table (i.e., IA, etc)
c     sequence = iso-electronic sequence (i.e., lithium, etc)
c     config   = iso-electronic sequence spectroscpic notation

      character*80      specID,ionID,shellID,sequence,group,config
      COMMON/chblk/     specID(Imax),ionID(Imax,Imax+1),
     &                  shellID(NSHmax),sequence(Imax),
     &                  group(Imax),config(Imax)


c     FILE PATH NAMES (file name generation, etc)

c     tabpath = directory path where atomic data is stored
c     UVBpath = directory path where UVB SED grids are stored
c     Sb99path = directory path where stellar pop SED grdis are stored

      character*250     tabpath,UVBpath,Sb99path
      COMMON/pathblk/   tabpath,UVBpath,Sb99path

