c     =====================================================================
c                                                                         
c               ART Version 3-1 : Read control parameters              
c                                                                         
c            by Alexei Khokhlov and Andrey Kravtsov (1996-1999)           
c                                                                         
c     =====================================================================
c
c     -------------------------
      subroutine ReadControl ()
c     -------------------------
c
#     include "a_def.h"
      INCLUDE 'a_tree.h'
      INCLUDE 'a_control.h'
c

      write (*,'(/1x, ''****************************************'',
     &           /1x, ''* ART N-body/Gas Version 3-1: 12-04-03 *'',
     &           /1x, ''****************************************'')')

      open ( unit = 20 , file = 'control.dat' , status = 'old' )
      rewind ( unit = 20 )
      read  ( unit = 20 , fmt = '(10x,    a256)')  jobname1
      read  ( unit = 20 , fmt = '(10x,    a256)')  path
      read  ( unit = 20 , fmt = '(70x,     l10)')  start
      read  ( unit = 20 , fmt = '(70x,     i10)')  irun
      read  ( unit = 20 , fmt = '(70x,     i10)')  nfsave
      read  ( unit = 20 , fmt = '(70x,     l10)')  lviewout
      read  ( unit = 20 , fmt = '(70x,     i10)')  ntc
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  gamma
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  cfl
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  timinc
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  rhohdr
      read  ( unit = 20 , fmt = '(70x,   e10.2)')  aHR

      radin(1) = gamma

#ifdef STARFORM
      open ( unit = 21 , file = 'sf.dat' , status = 'old' )
      rewind ( unit = 21 )
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  alpha_SF
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  eps_SF
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  dtmin_SF
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  dm_star_min
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  rho_SF
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  T_SF
      if     ( alpha_SF .le. 0.0 .or. alpha_SF .gt. 3.0 ) then 
        write(*,*) '* error in ReadControl: bad alpha_SF:', alpha_SF
        write(*,*) '* expected alpha_SF in the range (0.0,3.0)'
        stop
      elseif ( alpha_SF .lt. 0.99 .or. alpha_SF .gt. 2.01 ) then 
        write(*,*) 'warning: in SetUnits: suspicious alpha_SF =',
     &     alpha_SF
        write(*,*) 'hope this is what you wanted...'
      elseif ( eps_SF   .le. 0.0 ) then 
        write(*,*) '* error in ReadControl: bad eps_SF:', eps_SF
        write(*,*) '* eps_SF must be > 0'
        stop
      elseif ( dtmin_SF   .le. 0.0 .or. dtmin_SF .gt. 1.e11) then 
        write(*,*) '* error in ReadControl: bad dtmin_SF:', dtmin_SF
        write(*,*) '* expected dtmin_SF in the range (0.0,1e11)'
        stop
      elseif ( rho_SF   .lt. 0.0 ) then 
        write(*,*) '* error in ReadControl: bad rho_SF:', rho_SF
        write(*,*) '* rho_SF must be >= 0'
        stop
      elseif ( T_SF   .lt. 0.0 ) then 
        write(*,*) '* error in ReadControl: bad T_SF:', T_SF
        write(*,*) '* T_SF must be >= 0'
        stop
      endif
c
c.... if a_IMF > 0 -> power law IMF with the slope a_IMF
c     otherwise -> user specified IMF in f_IMF
c
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  a_IMF  
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  aM_stl
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  aM_stu
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  aM_SNII
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  aM_SNIa1
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  aM_SNIa2
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  ejM_SNIa
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  C_SNIa
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  t_SNIa
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  E_51
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  t_fb
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  T0_ml
      read  ( unit = 21 , fmt = '(70x,   e10.2)')  c0_ml

      t_SNIai = 1.d0 / t_SNIa
      if     ( aM_stl .ge. aM_stu ) then 
         call Open_ASCII_File ( 13 , errorfile , sequent )
        write(13,*) '1   : * error in ReadControl: bad aM_stl, aM_stu:',
     &      aM_stl, aM_stu
        write(13,*) '      * must be aM_stl < aM_stu '
        close ( 13 ) 
        stop
      elseif ( aM_SNII .ge. aM_stu ) then 
        call Open_ASCII_File ( 13 , errorfile , sequent ) 
        write(13,*)'1   : * error in ReadControl: bad aM_SNII, aM_stu:',
     &      aM_SNII, aM_stu
        write(13,*)'      * must be aM_SNII < aM_stu '
        close ( 13 ) 
        stop
      elseif ( aM_SNIa1 .ge. aM_SNIa2 ) then 
        call Open_ASCII_File ( 13 , errorfile , sequent ) 
        write(13,*)'1 : * error in ReadControl: bad aM_SNIa1,2:',
     &      aM_SNIa1, aM_SNIa2
        write(13,*)'      * must be aM_SNIa1 < aM_SNIa2 '
        close ( 13 ) 
        stop
      elseif ( aM_SNIa1 .ge. aM_stu ) then 
        call Open_ASCII_File ( 13 , errorfile , sequent ) 
        write(13,*)'1 : * error in ReadControl: bad aM_SNIa1, aM_stu:',
     &      aM_SNIa1, aM_stu
        write(13,*)'      * must be aM_SNIa1 < aM_stu '
        close ( 13 ) 
        stop
      elseif ( aM_SNIa2 .ge. aM_stu ) then 
        call Open_ASCII_File ( 13 , errorfile , sequent ) 
        write(13,*)'1 : * error in ReadControl: bad aM_SNIa2, aM_stu:',
     &      aM_SNIa2, aM_stu
        write(13,*)'      * must be aM_SNIa2 < aM_stu '
        close ( 13 ) 
        stop
      endif
#endif

!#ifdef REFINE
        read  ( unit = 20 , fmt = '(70x,   e10.2)') wsplit
        read  ( unit = 20 , fmt = '(70x,   e10.2)') wjoin
        read  ( unit = 20 , fmt = '(70x,   i10  )') MinL_Jeans 
        read  ( unit = 20 , fmt = '(70x,   i10  )') Ndiff 
        do iref = 1 , n_ref
          read(20,*) 
     &      irefdummy, w_ref(iref), floor_ref(iref),
     &      lmin_ref(iref), lmax_ref(iref)
        end do
        read (20,*) xrefmin, xrefmax
        read (20,*) yrefmin, yrefmax
        read (20,*) zrefmin, zrefmax
!#endif

      read(20,*) n_save
      if ( n_save .gt. 0 ) then 
        do is = 1 , n_save 
          read(20,*) asave(is)
        enddo
      endif
      

      close ( unit = 20 )

      write (*,'(/'' jobname1:'',/a256)') jobname1
      write (*,'( '' path    :'',/a256)') path

      write (*,'(/1x,''mcell   '', i10,10x,
     &               ''mind    '', i10,10x,
     &               ''NCPUS   '', i10)')   mcell, mind, NCPUS
      write (*,'( 1x,''ndim    '', i10,10x,
     &               ''nvar    '', i10,10x, 
     &               ''nchem   '', i10)')   ndim, nvar, nchem
      write (*,'( 1x,''start   '', l10,10x
     &               )')   start
      write (*,'( 1x,''ntc     '', i10,10x,
     &               ''cfl     '', 1pe10.2,10x,
     &               ''timinc  '', 1pe10.2)') ntc, cfl, timinc
      write (*,'( 1x,''rhogrv  '', 1pe10.2,10x,
     &               ''rhohdr  '', 1pe10.2,10x,
     &               ''enghdr  '', 1pe10.2)') rhogrv, rhohdr, enghdr
      write(*,*)
      write(*,*)
#ifdef STARFORM
      write (*,'( 1x,''alpha_SF  '', 1pe10.2,10x,
     &               ''eps_SF    '', 1pe10.2,10x,
     &               ''dtmin_SF  '', 1pe10.2)') alpha_SF,eps_SF,dtmin_SF
      write (*,'( 1x,''dm_st_min '', 1pe10.2,10x,
     &               ''rho_SF    '', 1pe10.2,10x,
     &               ''T_SF      '', 1pe10.2)') dm_star_min,rho_SF,T_SF
      write (*,'( 1x,''a_IMF     '', 1pe10.2,10x,
     &               ''M_l       '', 1pe10.2,10x,
     &               ''M_u       '', 1pe10.2)') a_IMF, aM_stl, aM_stu
      write (*,'( 1x,''M_SNII    '', 1pe10.2,10x,
     &               ''E_51      '', 1pe10.2,10x,
     &               ''t_fb      '', 1pe10.2)') aM_SNII, E_51, t_fb
      write (*,'( 1x,''M_SNIa1   '', 1pe10.2,10x,
     &               ''M_SNIa2   '', 1pe10.2,10x,
     &               ''C_SNIa    '', 1pe10.2)') 
     &                                      aM_SNIa1, aM_SNIa2, C_SNIa
      write (*,'( 1x,''t_SNIa    '', 1pe10.2,10x,
     &               ''ejM_SNIa  '', 1pe10.2)')
     &                                      t_SNIa, ejM_SNIa
#endif
#ifdef STELLARMASSLOSS
      write (*,'( 1x,''Stellar mass loss parameters :    '',
     &               ''c_ml      '', 1pe10.2,10x,
     &               ''T0_ml     '', 1pe10.2)') c0_ml, T0_ml 
#endif
!#ifdef REFINE
        write (*,'(/'' Refinement controls:'')')
        write (*,'(1x,''wsplit/wjoin            : '',
     &             1pe10.2,'' / '', 1pe10.2)') 
     &    wsplit, wjoin
        write (*,'(1x,''Ndiff   : '',
     &             i10)')
     &    Ndiff
        write (*,'(1x,''MinL_Jeans  : '',
     &             i10)')
     &    MinL_Jeans
        write (*,'(1x,''      iref'',
     &                ''     w_ref'',
     &                '' floor_ref'',
     &                ''  lmax_ref'')')
        do iref = 1 , n_ref
          write (*,'(1x,i10,1p2e10.2,2(1x,i5))') iref,
     &      w_ref(iref), floor_ref(iref), lmin_ref(iref), lmax_ref(iref)
        end do
        write(*,*) 'xrefmin,max =',xrefmin,xrefmax
        write(*,*) 'yrefmin,max =',yrefmin,yrefmax
        write(*,*) 'zrefmin,max =',zrefmin,zrefmax
!#endif

      write(*,'()')
      do i = 1 , n_save
        write(*,*) 'isave =',i,' asave =',asave(i)
      enddo

      write(*,'()')
       
      return
      end

