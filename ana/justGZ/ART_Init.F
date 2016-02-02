c     ===================================================================
c
c      ART Version 4: ART_Init.f - initialization of variables and arrays
c
c     ===================================================================
c
c     -------------------------
      subroutine Init_Arrays ()
c     -------------------------
c
      include 'a_tree.h'

CEVERINO08272009C$OMP PARALLEL SECTIONS
CEVERINO08272009C$OMP SECTION
C$OMP PARALLEL DO DEFAULT(SHARED) 
C$OMP+SCHEDULE(STATIC, 10000)
C$OMP+PRIVATE(j)
      do i = 1 , nctot
        iCL(i)    = nil 
        iOctCh(i) = nil 
        do j = 1 , nvar
          var(j,i) = zero
        enddo
        do j = 1 , nhvar
          hvar(j,i) = zero
        enddo
#ifdef RadPressure
          RadPre(i)=zero
#endif
       enddo

      do Level = MinLevel+1 , MaxLevel
        iHOLL(Level) = nil
        iNOLL(Level) = nil
      enddo

CEVERINO08272009C$OMP SECTION
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+SCHEDULE(STATIC, 10000)
C$OMP+PRIVATE(i)    
      do j = 1 , moct
        iOctLL1(j) = nil
        iOctLL2(j) = nil 
        iOctLv(j)  = iFreeLevel
CEVERINO08272009
        iOctPs(1,j)=nil
        iOctPs(2,j)=nil 
        iOctPs(3,j)=nil 
        do i=1,nneigh
           iOctNb(i,j)=nil
        enddo
        iOctPr(j)=nil
        iOctLv(j)=nil
        iOctCh(j)=nil

        iSelect(j)=nil
      enddo
CEVERINO08272009C$OMP END PARALLEL SECTIONS
      iOctFree     = 1
      nOct         = 0

C$OMP PARALLEL DO DEFAULT(SHARED)      
C$OMP+SCHEDULE(STATIC, 10000) 
      do j = 1 , npmax 
        iLL(j,1) = nil
        iLL(j,2) = nil 
      enddo

      return
      end
c
c     ---------------------
      subroutine SF_Init ()
c     ---------------------
c
c     initialize starformation/feedback/enrichment related constants
c
#     include "a_def.h"
      include 'a_numbers.h'
      include 'a_setup.h'
      include 'a_control.h'

      real*8 dummy1, dummy2, dummy3, dummy4, dummy5, dumst
      real*8 fp(10), Aprime
      real*8 INTEGRATE, f_IMF, fm_IMF, f_IMF_plaw, fm_IMF_plaw
      real*8 fej_IMF, fej_IMF_plaw, f_SNIa_func
      real*8 aN_SNII, An_IMF
      external f_IMF, fm_IMF, f_IMF_plaw, fm_IMF_plaw
      external fej_IMF, fej_IMF_plaw, f_SNIa_func
c
      open (iOIO, file = TRIM(iofile), position='append')
c
#ifdef ENRICH
      if ( nhvar .lt. 9 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*) 
     &     '1  : SF_Init : ENRICH is defined in a_def.h but'
          write(13,*)' nhvar =', nhvar,' < 9 in a_setup.h'
          write(13,*)'stopping...'
          close ( 13 ) 
          stop          
      endif
#endif      
#ifdef ENRICH_SNIa
#ifdef ENRICH
      if ( nhvar .lt. 10 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*) 
     &     '1  : SF_Init: ENRICH and ENRICH_SNIa are defined in a_def.h'
          write(13,*)' but nhvar =', nhvar,' < 10 in a_setup.h'
          write(13,*)'stopping...'
          close ( 13 ) 
          stop          
      endif
#endif
#ifndef ENRICH
      if ( nhvar .lt. 9 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*) 
     &     '1  : SF_Init : ENRICH_SNIa is defined in a_def.h'
          write(13,*)' but nhvar =', nhvar,' < 9 in a_setup.h'
          write(13,*)'stopping...'
          close ( 13 ) 
          stop          
      endif
#endif
#endif      
c
c.... below t0 must be in yrs, rho0 in Msnu / Mpc^-3 (no h!) 
c
c
c.... starformation efficiency factor c_SFR
c
      C_SFR = eps_SF * 2.5*10.0**(6.0 - 16.0*alpha_SF) * 
     &        t0 * rho0**(alpha_SF-1.0)

      E_0 = 1.38d58 * r0**5 * Om0**2 / hubble  ! [ergs]
      write(iOIO,*) 'E_0, r0, Om0, aM0, sun_mass =',
     &     E_0, r0, Om0, aM0, sun_mass  

      IF ( a_IMF .gt. 0.0 ) THEN  ! if power-law
        if ( a_IMF .eq. 2.0 ) then 
          Aprime = 1.0 / log(aM_stu/aM_stl)
        else
          Aprime = (2.0 - a_IMF) / 
     &             (aM_stu**(2.0-a_IMF) - aM_stl**(2.0-a_IMF)) 
        endif
        An_IMF = Aprime ! IMF normalization
        if ( a_IMF .eq. 1.0 ) then 
          Aprime = Aprime * log(aM_stu/aM_SNII)
        else
          Aprime = Aprime /
     &             (1.0 - a_IMF) * 
     &             (aM_stu**(1.0-a_IMF) - aM_SNII**(1.0-a_IMF)) 
        endif
#ifdef ENRICH
c
c....   initialize fraction of stellar mass returned as metals
c
        if ( a_IMF .eq. 2.0 ) then 
          dummy1 = 1.d0 / log(aM_stu/aM_stl)
        else
          dummy1 = (2.d0 - a_IMF) / 
     &             (aM_stu**(2.0-a_IMF) - aM_stl**(2.0-a_IMF)) 
        endif
        fp(1) = a_IMF
        dumst = 1.d-1 * (aM_stu - aM_SNII)
        dummy2 = 
     &    INTEGRATE(fej_IMF_plaw, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        fmass_met = dummy2 * dummy1 
#endif
        fp(1) = a_IMF
        dummy2 = 
     &    INTEGRATE(fm_IMF_plaw, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        dummy1 = 
     &    INTEGRATE(f_IMF_plaw, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)

        aMSN_ave = dummy2 / dummy1
        aN_SNII = dummy1 * An_IMF
        
        dumst = 1.d-1 * (aM_stu - aM_stl)
        dummy2 = 
     &    INTEGRATE(fm_IMF_plaw, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        dummy1 = 
     &    INTEGRATE(fm_IMF_plaw, fp, 1, aM_stl, aM_stu, dumst, 1.d-9)
        fmass_SN = dummy2 / dummy1 

#ifdef FEEDBACK_SNIa         
        dummy1 = 
     &   INTEGRATE(fm_IMF_plaw, fp, 2, aM_SNIa1, aM_SNIa2, dumst, 1.d-9)
        dummy2 = 
     &   INTEGRATE(f_IMF_plaw, fp, 2, aM_SNIa1, aM_SNIa2, dumst, 1.d-9)
        dummy3 = 
     &   INTEGRATE (f_SNIa_func, fp, 2, 1.d-2, 1.d3 , 1.d-1 , 1.d-9)
c
c....   to get SNIa rate in SN/yr: R_Ia = Mstar * 1.d-9 * RIaf * f_SNIa
c
        RIaf = 1.d-9 * t0 * aM0 * C_SNIa * dummy2 / dummy1 / dummy3 /
     &       t_SNIa
#endif
        
      ELSE ! Miller-Scalo IMF
        dumst = 1.d-1 * (aM_stu - aM_stl)
        dummy1 = 
     &    INTEGRATE ( fm_IMF, fp, 1, aM_stl, aM_stu, dumst, 1.d-9 )
        An_IMF = 1.d0 / dummy1

        dumst  = 1.d-1 * (aM_stu - aM_SNII)
        dummy2 = 
     &    INTEGRATE ( f_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9 )
        if ( dummy1 .le. 0.0 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*) 
     &     '1  : SF_Init : something is wrong dummy1 <= 0 :', 
     &     dummy1
          close ( 13 ) 
          stop          
        endif
c
        Aprime = dummy2 / dummy1
c
        if ( Aprime .le. 0.0 ) then 
          call Open_ASCII_File ( 13 , errorfile , sequent ) 
          write(13,*)' 1  : SF_Init : something is wrong Aprime <= 0 :',
     &      Aprime
          close ( 13 ) 
          stop          
        endif 
#ifdef ENRICH
c
c....   initialize fraction of stellar mass returned as metals
c
        
        dumst = 1.d-1 * (aM_stu - aM_SNII)
        dummy2 = 
     &    INTEGRATE (fej_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)        
        fmass_met = dummy2 / dummy1         
#endif
        dummy2 = 
     &    INTEGRATE(fm_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        dummy1 = 
     &    INTEGRATE(f_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)

        aMSN_ave = dummy2 / dummy1
        aN_SNII = dummy1 * An_IMF
 
        dumst = 1.d-1 * (aM_stu - aM_stl)
        dummy2 = 
     &    INTEGRATE(fm_IMF, fp, 1, aM_SNII, aM_stu, dumst, 1.d-9)
        dummy1 = 
     &    INTEGRATE(fm_IMF, fp, 1, aM_stl, aM_stu, dumst, 1.d-9)
        fmass_SN = dummy2 / dummy1 
#ifdef FEEDBACK_SNIa
        dumst = 1.d-1 * ( aM_SNIa2 - aM_SNIa1 )
        dummy1 = 
     &    INTEGRATE ( fm_IMF, fp, 2, aM_SNIa1, aM_SNIa2, dumst, 1.d-9 )
        dummy2 = 
     &    INTEGRATE ( f_IMF, fp, 2, aM_SNIa1, aM_SNIa2, dumst, 1.d-9 )
        dummy3 = 
     &    INTEGRATE ( f_SNIa_func, fp, 2, 1.d-3, 1.d3 , 1.d-1 , 1.d-9)
c
c....   to get SNIa rate R_Ia =  Mstar * RIaf * f_SNIa(dt)
c


        RIaf = 1.d-9 * t0 * aM0 * 
     &         C_SNIa * dummy2 / dummy1 / dummy3 / t_SNIa
#endif
      ENDIF

#ifdef ENRICH
      write(iOIO,*) 'fmass_met =',fmass_met
#endif
      write(iOIO,*) 'fraction of mass in SNII:', fmass_SN
      write(iOIO,*) 'average M_SN=', aMSN_ave
      write(iOIO,*) 'N of SNII:', aN_SNII
      write(iOIO,*) 'N of SNIa:', C_SNIa  * dummy2/dummy1

      C_fb = 1.d51 * E_51 * Aprime * aM0 / E_0       
      C_fbIa = 1.d51 * E_51 / E_0       

      rho_SF_fact = rho_SF * 8.9d4 / hubble**2 / Om0 / (1.d0 - Y_p)
      write(iOIO,*) 'C_fb, C_fbIa, E_51, Aprime, aM0, E_0, rho_SF =',
     &      C_fb, C_fbIa, E_51, Aprime, aM0, E_0, rho_SF_fact
c
      close (iOIO)
c
      return
      end


c     ----------------------------------------
      real*8 function f_SNIa_func ( xd , fp , np )
c     ----------------------------------------
c
c     function for integrating S R_SNIa(t) dt
c
c     xd = t_SNIa / t
c
      integer np
      real*8 xd, fp(np)

      f_SNIa_func = max ( dexp(-xd*xd) / dsqrt(xd) , 1.d-20 )

      return
      end

c     ----------------------------------------
      real*8 function f_SNIa ( xd , fp , np )
c     ----------------------------------------
c
c     function describing phenomenological shape of the
c     R_SNIa(t) function
c
c     x = t_SNIa / t
c
      integer np
      real*8 xd, fp(np)

      f_SNIa = max ( dexp(-xd*xd) * xd**1.5 , 1.d-20 )

      return
      end

c     -----------------------------------------------
      real*8 function f_IMF_plaw ( amstar , fp , np )
c     -----------------------------------------------
c
c     unnormalized stellar IMF approximated by a single power law
c
c     input: real*8 amstar - star mass in M_sun
c            fp(1) - power index (e.g., 2.35 for the Salpeter IMF)

      integer np
      real*8 amstar, fp(np)
c
      f_IMF_plaw = amstar**(-fp(1))

      return
      end

c     ------------------------------------------------
      real*8 function fm_IMF_plaw ( amstar , fp , np )
c     ------------------------------------------------
c
c     unnormalized stellar IMF approximated by a single power law
c         multiplied by mass
c
c     input: real*8 amstar - star mass in M_sun
c            fp(1) - power index (e.g., 2.35 for the Salpeter IMF)

      integer np
      real*8 amstar, fp(np)
c
      fm_IMF_plaw = amstar**(1.d0 - fp(1))

      return
      end

c     ------------------------------------------------
      real*8 function fej_IMF_plaw ( amstar , fp , np )
c     ------------------------------------------------
c
c     unnormalized stellar IMF approximated by a single power law
c         multiplied by mass and by fraction of mass ejected as metals
c
c     input: real*8 amstar - star mass in M_sun
c            fp(1) - power index (e.g., 2.35 for the Salpeter IMF)

      integer np
      real*8 amstar, fp(np)
      real*8 fmet_ej
c
      fej_IMF_plaw = amstar**(1.d0 - fp(1)) * fmet_ej ( amstar )

c      write(*,*) amstar, fmet_ej ( amstar )

      return
      end

c     ------------------------------------------
      real*8 function f_IMF ( amstar , fp , np )
c     ------------------------------------------
c
c     unnormalized stellar IMF
c
c     input: real*8 amstar - star mass in M_sun
c            fp(np) - dummy array, is not used here

      integer np
      real*8 amstar, fp(np)
      real*8 C_1, C_2
      parameter ( C_1 = 1.09 , C_2 = -1.02 )   ! parameters of ther M-S IMF
      real*8 am0_Ch, beta_Ch, alpha_Ch
      parameter ( am0_Ch = 716.4 , beta_Ch = 0.25, alpha_Ch = -3.3 )

c========= STG: switched from M-S to Chabrier:
c
c.... Miller-Scalo (1979, ApJS 41, 513, eq. 30, Table 7) IMF
c
c      f_IMF = dexp ( -C_1 * (dlog10(amstar) - C_2)**2 ) / amstar
c
c.... Chabrier, G. (2001, ApJ 554, 1274)
c
      f_IMF = dexp(-(am0_Ch/amstar)**beta_Ch) * amstar**alpha_Ch
c===========

      return
      end

c     ------------------------------------------
      real*8 function fm_IMF ( amstar , fp , np )
c     ------------------------------------------
c
c     unnormalized stellar IMF * m_star
c
c     input: real*8 amstar - star mass in M_sun
c            fp(np) - dummy array, is not used here

      integer np
      real*8 amstar , fp(np), f_IMF

      fm_IMF = amstar * f_IMF ( amstar , fp , np )

      return
      end

c     ------------------------------------------
      real*8 function fej_IMF ( amstar , fp , np )
c     ------------------------------------------
c
c     unnormalized stellar IMF * m_star * fraction of mass ejected as metals
c
c     input: real*8 amstar - star mass in M_sun
c            fp(np) - dummy array, is not used here

      integer np
      real*8 amstar , fp(np), f_IMF, fmet_ej

      fej_IMF = amstar * f_IMF ( amstar , fp , np ) * fmet_ej ( amstar )

      return
      end

c     ----------------------------------
      real*8 function fmet_ej ( amstar )
c     ----------------------------------
c
c     for amstar >= M_SN, returns fraction of mass ejected as metals into ISM
c
c     current approximation is based on the results of Woosley & Weaver 1995
c     and ignores metallicity dependence of the yields
c
c     input : real*8 amstar in Msun
c
      real*8 amstar

      fmet_ej = min ( 2.0d-1, max(1.0d-2*amstar - 6.d-2,1.d-20) )

      return
      end


