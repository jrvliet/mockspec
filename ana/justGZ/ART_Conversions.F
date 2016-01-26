c     ====================================================================
c                                                                         
c        HART Version 1-0 : various routine converting to and from        
c                           internal code units
c                                                                         
c                      by Andrey Kravtsov (1998-1999)           
c                                                                            
c     =====================================================================
c     
c     ----------------------
      subroutine SetUnits ()
c     ----------------------
c
c     set units for physical variables
c
c
#     include "a_def.h"
      include 'a_numbers.h'
      include 'a_setup.h'
      include 'a_control.h'
      real wmu53
c
      if ( boxh .le. 0. ) then 
        write(*,*) '* SetUnits : something is wrong: boxh <= 0 :',boxh
        stop
      endif
      wmu53 = wmu**(5./3.)
c
c.... comoving size of the zeroth-level cell in /h Mpc
c
      r0 = boxh / ng  
c
      if ( Om0 .le. 0. ) then 
        write(*,*) '* SetUnits : something is wrong: Om0 <= 0 :',Om0
        stop
      endif
      if ( hubble .le. 0. ) then 
        write(*,*)'* SetUnits : something is wrong: hubble <= 0:',hubble
        stop
      endif
      if ( hubble .gt. 1. ) then 
        write(*,*)'* warning : SetUnits : suspicious hubble =',hubble
      endif
c        
      dummy = 100.0 * hubble * sqrt(Om0) 
c
c.... time unit in yrs
c 
      t0 = 2.0 / dummy * 3.0856e19 / 3.15e7
c
c.... v0 - velocity units in km/s
c
      v0 = 50. * r0 * sqrt(Om0)
c
c.... rho0 = 3H0^2 * Om0 / (8*pi*G) - unit of density in Msun/Mpc^3
c
      rho0 = 2.776e11 * hubble**2 * Om0
c
c.... P0 = rho0 * v0**2 - unit of pressure in g/cm/s^2
c
      P0 = 4.697e-16 * Om0**2 * r0**2 * hubble**2
c
c.... T_0 = unit of temperature in K and in keV)
c
      T_0 = 3.03e5 * r0**2 * wmu * Om0 ! [K]
c      T_0 = 2.61155 * r0**2 * wmu * Om0 ! [keV]
c
c.... S_0 = unit of entropy in keV * cm^2
c
      S_0 = 52.077 * wmu53 * hubble**(-4./3.)*Om0**(1./3.)*r0**2
c
c.... mass conversion (Mbox = rho0 * Lbox^3, Mbox_code = Ng^3
c     for non-cosmological run aM0 must be defined during initialization
c     [aM0] = [Msun]
      aM0 = rho0 * (boxh/hubble)**3 / ncell0
c
c.... constant L0 for the Sutherland & Dopita cooling function
c
c.... old
c      AL_SD = 2.0731d+22 * (4.d0 - 3.d0*Y_p) * (2.d0 - Y_p) * hubble /
c     &        sqrt(Om0) / r0**2  ! in ergs/s/cm^3
c.... new
c
       AL_SD = 1.6625d0 * (1.d0 - Y_p)**2 * hubble /
     &         sqrt(Om0) / r0**2  ! in ergs/s/cm^3
c
c.... constant for Compton cooling
c
      AL_Comp = 4.55d-4 * (4.d0 - 2.d0 * Y_p) / (8.d0 - 5.d0 * Y_p) * 
     &            T_CMB0**4 / hubble / sqrt(Om0)
c
#ifdef STARFORM 
      call SF_Init () 
#endif
      open (iOIO, file = TRIM(iofile), position='append')
      write (iOIO,*) 'set units to :'
      write (iOIO,'(1x,''r0 [h^{-1} Mpc]    = '',1pe12.4)') r0
      write (iOIO,'(1x,''t0 [yrs]           = '',1pe12.4)') t0
      write (iOIO,'(1x,''rho0 [M_sun/Mpc^3] = '',1pe12.4)') rho0
      write (iOIO,'(1x,''v0 [km/s]          = '',1pe12.4)') v0
      write (iOIO,'(1x,''T0 [K]             = '',1pe12.4)') T_0
      write (iOIO,'(1x,''P0 [g/cm/s^2]      = '',1pe12.4)') P0
      write (iOIO,'(1x,''S0 [keV cm^2]      = '',1pe12.4)') S_0
#ifdef STARFORM
      write (iOIO,'(1x,''E0 [ergs]          = '',1pe12.4)') E_0
      write (iOIO,'(1x,''C_fb               = '',1pe12.4)') C_fb
      write (iOIO,'(1x,''C_sfr              = '',1pe12.4)') C_SFR      
      write (iOIO,'(1x,''rho_SF             = '',1pe12.4)') rho_SF
#endif
      close (iOIO)
      return
      end
c
c     --------------------------
      real*8 function a2b ( at )
c     --------------------------
c
c     translates expansion factor at into hydro time variable b
c
      real*8 at
      real*8 fp(2), atst
      real*8 INTEGRATE , f_a2b
      external INTEGRATE , f_a2b

      include 'a_numbers.h'
      include 'a_setup.h'
      include 'a_control.h'
c      
c
      if ( (at .lt. zero) ) then 
        write(*,*) 'a2b : something is wrong : at = ', at
        stop
      endif
      IF ( (Om0 .eq. one) .and. (Oml0 .eq. zero) ) THEN 
        a2b = (1.d0 - 1.d0/sqrt(at))
      ELSE
        if ( (Om0 .eq. zero) .and. ( Oml0 .eq. zero ) ) then 
          write(*,*) 'a2b : something is wrong : Om0 = Oml0 = 0.0'
          stop
        endif
        fp(1) = Om0
        fp(2) = Oml0
        atst  = 1.d-1 * (1.d0 - at)
        a2b = INTEGRATE( f_a2b , fp , 2 , 1.d0 , at , atst , 1.d-9 )
      ENDIF       
c
      return
      end
c
c     -------------------------------------
      real*8 function f_a2b ( x , fp , np )
c     -------------------------------------
c
c     input : x     - is expansion factor variable
c             fp(1) - Om0  = present-day matter density 
c             fp(2) - Oml0 = present-day vacuum contribution
c
      integer np 
      real*8 x , fp(np), d1
c
      d1 = x**3
      f_a2b = .5d0 * sqrt(fp(1)) / d1 
     &             / sqrt(fp(1)/d1 + fp(2) + (1.d0-fp(1)-fp(2))/x**2)
c
      return
      end
c
c     --------------------------
      real*8 function b2a ( bt )
c     --------------------------
c
      real*8 bt
      real*8 f_b2a, fp(1)
      real zbrent 
      external zbrent, f_b2a
c
      include 'a_numbers.h'
      include 'a_setup.h'
      include 'a_control.h'
c
      IF ( (Om0 .eq. one) .and. (Oml0 .eq. zero) ) THEN 
        b2a = (1.d0 / (1.d0 - bt))**2
      ELSE
        if ( (bt .lt. -190.d0) .or. (bt .gt. 1.d-1) ) then 
          write(*,*) ' b2a : something is wrong : bt =',bt
          stop
        endif
        if ( (Om0 .eq. zero) .and. ( Oml0 .eq. zero ) ) then 
          write(*,*) 'b2a : something is wrong : Om0 = Oml0 = 0.0'
          stop
        endif
        fp(1) = bt
        b2a = zbrent ( f_b2a , fp , 1 , 1.d-4 , 1.1d0 , 1.d-9 )
      ENDIF
c
      return
      end
c
c     --------------------------------------
      real*8 function f_b2a ( at , fp , np )
c     --------------------------------------
c
c     at - expansion factor
c     fp(1) = b
c
      integer np
      real*8 at , fp(np)
      real*8 a2b
      external a2b
c
      f_b2a = a2b ( at ) - fp(1)
c      
      return
      end
c
c     ---------------------------
      real*8 function age ( td )
c     ---------------------------
c
c     returns age of the Universe at td (time in code units)
c
c     uses Oleg's formula for flat LCDM 
c
      include 'a_tree.h'
      include 'a_control.h'
c
      real*8 td, b2a

      as = b2a ( td ) ! convert code time to expansion factor
      ff = Om0/(1.-Om0)/as**3 
c
c  .... calculate age of the universe in Gyrs 
c
      age = 9.779/hubble * 2./3./
     &      sqrt(1.-Om0)*log((1.+sqrt(1.+ff))/sqrt(ff))

      return
      end

