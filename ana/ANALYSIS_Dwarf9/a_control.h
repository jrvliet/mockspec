c     --------------------------------
c     Control parameters and variables
c     --------------------------------
c
c======STG:
      integer iExitCode2
      common / RUN2 / iExitCode2
c      logical Continue_Run
c=============
c======STG:
c      integer nstep              
c      common / SEBASTIAN / nstep 
c=============
c
c
      parameter ( nextra = 2 )  ! number of additional parameters 
      real extra(nextra)
      character*256 lextra(nextra)
      common / EXTRA01 / lextra
      common / EXTRA02 / extra     

      real boxh, Om0, Oml0, Omb0, hubble, aexpn, ainit, gamma 
      common / RUNPARAM / boxh, Om0, Oml0, Omb0, hubble,
     &                    aexpn, ainit, gamma

      real*8 r0, rho0, v0, t0, T_0, P0, S_0, aM0, E_0
      common / UNITS01 / r0, 
     &                   rho0, 
     &                   v0, 
     &                   t0, 
     &                   T_0, 
     &                   P0, 
     &                   S_0, 
     &                   aM0, 
     &                   E_0
      real*8 AL_SD, AL_Comp
      common / COOL01 / AL_SD, AL_Comp
c
c     ....... starformation parameters ........
c
      real*8 alpha_SF, C_SFR, eps_SF, dtmin_SF, dm_star_min, 
     &       rho_SF, rho_SF_fact, T_SF, a_IMF, aM_stl, aM_stu, aM_SNII, 
     &       aM_SNIa1,aM_SNIa2, t_SNIa, t_SNIai, C_SNIa, RIaf, ejM_SNIa,
     &       E_51, t_fb, C_fb, C_fbIa, fmass_met, c0_ml, T0_ml 
      common / SF01 / alpha_SF, C_SFR, eps_SF, dtmin_SF, dm_star_min, 
     &       rho_SF, rho_SF_fact, T_SF, a_IMF, aM_stl, aM_stu, aM_SNII, 
     &       aM_SNIa1,aM_SNIa2, t_SNIa, t_SNIai, C_SNIa, RIaf, ejM_SNIa,
     &       E_51, t_fb, C_fb, C_fbIa, c0_ml, T0_ml 
      common / SF03 / fmass_met
c
c     indices of hvar for metal products of SN type II & Ia
c
      integer izII, izIa  
      parameter ( izII = 9 , izIa = 10 ) 
c      
c     ............numerical controls...........
c
      common /NUMC01/ niter                      ! # of relax. iterations
c
c     ........interpolation coefficients.......
c
      common /INTER1/  wa, wbcd                  ! prolongation weights
c      
c     ....run control parameters & variables...
c
      real AEXP0,AMPLT,ASTEP,PARTW,
     &               TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,ap1,ap0,
     &               Wp5,extras(100)
c
      common / RUN / AEXP0,AMPLT,ASTEP,PARTW,
     &               TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,ap1,ap0,Wp5,
     &               NROWC,NGRIDC,nspecies,Nseed,istep,istep2,extras
c
c     ..................time...................

      real*8 t, dtime0, dtmin
      integer istep2
      real*8 tl(MinLevel:MaxLevel), tlold(MinLevel:MaxLevel)
      real*8 dtl(MinLevel:MaxLevel), dtlold(MinLevel:MaxLevel)
      real*8 aexp(MinLevel:MaxLevel), aexpold(MinLevel:MaxLevel)
      integer iTimeBin(MinLevel:MaxLevel)	  
      integer iSO(MinLevel:MaxLevel)  ! sweep order
c
      common / TIME01 /  t       ,   ! current time of the simulation
     &                   dtime0  ,   ! current time step
     &                   tl      ,   ! current  time moment for level L
     &                   tlold   ,   ! previous time moment for level L
     &                   dtl     ,   ! = dtime0/2**iTimeBin
     &                   dtlold  ,   ! previous time step for level L
     &                   dtmin   ,   ! = dtime(MaxLevel)
     &                   aexp    ,   ! = expansion factor (Level)
     &                   aexpold     ! = previous expansion factor (Level)
      common / TIME03 /  iTimeBin, iSO
c
c
c     .........refinement controls...........
c
      common /RFMC01/ nsplit(MinLevel:MaxLevel)   ! number  thresholds 
      common /RFMC02/ trho  (MinLevel:MaxLevel),  ! density thresholds 
     &                tbmass(MinLevel:MaxLevel),
     &                tphi  (MinLevel:MaxLevel)   ! potential threshold
c
c              trho = wpar * nsplit
c
c     .........refinement weights............
c
      parameter ( n_ref = 9 ) 
      integer lmin_ref(n_ref), lmax_ref(n_ref) ! min & max ref. level for each refinement index
      common / REFW /  wsplit , wjoin , w_ref(n_ref)
c
c     .....refinement density thresholds.....
c
      common / REFD01 /  floor_ref(n_ref)
      common / REFD02 /  lmin_ref, lmax_ref
c
c.... min/max of x, y, z defining rect. region where refinement is allowed
c  
      common / REFREGION / xrefmin, xrefmax, 
     &                     yrefmin, yrefmax, 
     &                     zrefmin, zrefmax
c
c
c     ..........tree related controls..........
c
      common /CELLS1/  ncell,                    ! number of cells in use
     &                 MaxLevelNow               ! current maximum level
      common /CUBE01/  xmin , ymin , zmin,       ! coordinates of cube corners
     &                 xmax , ymax , zmax
c
      character*45      HEADER
      common / HEADDR/  HEADER
c
c     .................timing..................
c
      parameter ( ntiming = 20 ) 
      common / TIMECPU / CPU(ntiming)
c
      parameter ( iostep = 1  )  ! output step
      parameter ( nsave  = 200 )
c
      dimension   isave (nsave)
      dimension   asave (nsave)                  ! time moments to save
      character*5 aname (nsave)
      common / SAVE01 / n_save, isave, asave, aname
c
c.... general controls
c
      logical start, 
     &        lviewout
c
      character*256 path
      character*256 jobname1
c
      real*8 enghdr
      common /CNTRL/   start,        ! T/F = start/continue
     &                 irun,         ! input/output routines
     &                 mstep,        ! make this number of time steps
     &                 nprint,       ! print step
     &                 nfsave,       ! dump step
     &                 nplot,        ! plot step
     &                 ntc,          ! increase dtime only after ntc steps
     &                 tremin,       ! exit if time left <= tremin
     &                 cfl,          ! Courant Number (<=1)
     &                 timinc,       ! time step increase factor
     &                 lhydro,       ! T/F = Euler flux on/off
     &                 lapidus,      ! T/F = Lapidus artif. diffusion on/off
     &                 lmassdiff,    ! T/F = artificial mass diffusion on/off
     &                 lvisc,        ! T/F = viscosity on/off
     &                 lcooling,     ! T/F = cooling on/off
     &                 lcond,        ! T/F = thermal diffusion on/off
     &                 ldiff,        ! T/F = species diffusion on/off
     &                 lchem,        ! T/F = chemistry on/off
     &                 lgrav,        ! T/F = gravity on/off
     &                 lnbody,       ! T/F = particles on/off
     &                 lcosmology,   ! T/F = cosmological run/non-cosmological
     &                 lviewout,     ! T/F = output for view? 
     &                 rhohdr,       ! density floor for the Magic subroutine
     &                 lref,         ! T/F = refinement on/off
     &                 Ndiff,        ! number of diffusion interations
     &                 MinL_Jeans,    ! min. Level for Jeans limit on int. energy
     &                 aHR           ! a value at which another level (L=8) is allowed at high redshift (CEVERINO05222008)

c
      common /CNTRL1/ path, jobname1
c
      common / SEQUENCE / imoviestep
c
      parameter ( n_nodes_max = 500 )
      character*256 directory(n_nodes_max)

      integer i_node, n_nodes
      common / NODES1 / i_node, n_nodes, directory

      integer iOWork, iOErr, iODebug, iOut
      parameter ( iOut = 6 , iOWork = 12 , iOErr = 13 , iODebug = 14 ) 

      character*256 workfile, errorfile, processfile
      character*256 cflfile, runfile, timingfile
      character*256 debugfile, starfile, agnfile


      parameter ( runfile = 'run.log ' )
      parameter ( workfile = 'work.log ' )
      parameter ( errorfile = 'error.log ' )
      parameter ( processfile = 'process.log ' )
      parameter ( timingfile = 'timing.log ' )
      parameter ( cflfile = 'cfl.log ' )
      parameter ( starfile = 'sf.log ' )
      parameter ( debugfile = 'debug.log ' )
      parameter ( agnfile = 'agn.log' )

      character*256 append, rewind, sequent
      parameter ( append = 'append ' )
      parameter ( rewind = 'rewind ' )
      parameter ( sequent = 'sequential ' ) 


