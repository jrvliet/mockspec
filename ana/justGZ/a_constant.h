c     ------------
c     a_constant.h       Universal constants                Units
c     ------------

c      implicit real*8 (a-f,o-z)

      parameter ( pi         = 3.14159265358978245      )
      parameter ( pi2        = 2.0 * pi                 )
      parameter ( pi4        = 4.0 * pi                 ) 
      parameter ( pi43       = pi4 / 3.0                )
      parameter ( clight     = 2.99792458e10            ) ! [cm/s]
      parameter ( planck     = 6.62618e-27              ) ! [erg s]
      parameter ( planck_    = planck / (2.0*pi)        ) 
      parameter ( boltz      = 1.380662e-16             ) ! [erg/K]
      parameter ( bolr       = 1.0 / boltz              )
      parameter ( avogadro   = 6.022045e23              ) ! [1/mole]
      parameter ( ravogadro  = 1. / avogadro            ) 
      parameter ( e_charge   = 4.80286e-10              ) ! [esu]
      parameter ( emev       = 1.6021892e-6             ) ! ?
      parameter ( e_mass     = 9.10953e-28              ) ! [g]
      parameter ( p_mass     = 1.67262e-24              ) ! [g]
      parameter ( grav_c     = 6.67259e-8               ) ! [dyn cm^2 / g^2]
      parameter ( e_energy   = e_mass * clight * clight )
      parameter ( emc2       = 2.0 * e_energy           ) 
      parameter ( sun_mass   = 1.991e33                 ) ! [g]
      parameter ( sun_radius = 6.960e10                 ) ! [cm]

      real*8 mpc, gyr
      parameter ( mpc        = 3.0856d19                ) ! km
      parameter ( gyr        = 3.1558d16                ) ! s
c
      real wmu
      parameter ( Y_p = 0.245 ) ! He mass fraction
      parameter ( wmu = 4.0 / (8.0 - 5.0 * Y_p) ) ! mol weight
      parameter ( T_CMB0 = 2.726 ) 
      parameter ( T_min = 300.0 ) ! T floor in K
