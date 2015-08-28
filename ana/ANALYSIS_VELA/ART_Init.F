c     ===================================================================
c
c      ART Version 3: ART_Init.f - initialization of variables and arrays
c
c     ===================================================================

c     -----------------------------
      subroutine Init_Parameters ()
c     -----------------------------
c
#     include "a_def.h"
      include 'a_tree.h'
      include 'a_control.h'

      real    wspecies(nspec)
      integer lspecies(nspec)
      equivalence (wspecies(1),extras(1)), (lspecies(1),extras(11))
c
#ifdef ENRICH
      if ( nhvar .le. 8 ) then 
        call Open_ASCII_File ( 13 , errorfile , sequent ) 
        write(13,*) '1: error in Init_Parameters :'
        write(13,*) 'ENRICH is defined in a_def.h but no storage is'
        write(13,*) 'dedicated in a_setup.h: need nhvar >=9'
        write(13,*) 'but nhvar=',nhvar        
        close ( 13 ) 
        stop
      endif
#endif      

c.... read in cosmology parameres if run is cosmological
c
c.... define cube corners

      xmin = 1.0
      ymin = 1.0 
      zmin = 1.0
      xmax = float(ng) + 1.0
      ymax = float(ng) + 1.0
      zmax = float(ng) + 1.0

      istep2 = 0

c      aexpn = .1
c      astep = 0.02

c.... initialize parameters for pyramidal interpolation
c      wa   = sqrt(11.)/(sqrt(11.) + 3.*sqrt(3.))   ! 0.38
c      wbcd = sqrt(3.)/(sqrt(11.) + 3.*sqrt(3.))

      wa   = 0.25
      wbcd = (1. - wa)/3

      MaxLevelNow = nil
      

c.... initialize level-dependent density threshold

c      nsplit(MinLevel) = 5
c      trho(MinLevel) = wpar(1) * nsplit(MinLevel) - 1.0  
c      do Level = MinLevel + 1 , MaxLevel
c        if ( Level .le. 1 ) then 
c          nsplit(Level) = 5
c        else
c          nsplit(Level) = 5
c        endif
c        trho(Level) = wpar(1) * nsplit(Level)
c      enddo

c.... initialize cell sizes

      do Level = MinLevel , MaxLevel
        CellSize(Level)  = 1.0 / 2**Level        
        CellSizei(Level) = 2**Level        
        CellVol (Level)  = 1.0 / 2**(3*Level)
        CellVoli(Level)  = 2.0**(3*Level)
      enddo

      do Level = MinLevel , MaxLevel 
        iTimeBin(Level) = Level 
      enddo
c
      return
      end

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

c     -----------------------
      subroutine Init_Tree ()  
c     -----------------------
c
c     initializes tree arrays & variables
c 
      include 'a_tree.h'
c
C$OMP SINGLE
      do i = 1 , moct
        if ( iOctLv(i) .eq. iFreeLevel .and. i .lt. moct ) then 
          inext = i+1
          do while (iOctLv(inext).ne.iFreeLevel .and. inext.lt.moct) 
            inext = inext + 1
          enddo
          if ( inext .lt. moct ) then 
             iOctPr(i) = inext  ! establish linked list of Octs #>0
          else
            if ( iOctLv(inext) .eq. iFreeLevel ) then 
              iOctPr(i) = inext
              iOctPr(inext) = 0
            else
              iOctPr(i) = 0 
            endif
          endif
        else
          if ( i .eq. moct ) then 
            iOctPr(i) = 0
          endif
        endif
      end do
C$OMP END SINGLE
c
      return
      end
c

c     ----------------------------
      subroutine Init_Cosmology ()
c     ----------------------------
c
c     Read in primary and initialize secondary cosmological parameters
c
      include 'a_tree.h'
      include 'a_control.h' 
c
      open ( 30 , file = 'cosmology.dat' ) 
c
      read  ( unit = 30 , fmt = '(70x,   e10.2)')  Om0
      read  ( unit = 30 , fmt = '(70x,   e10.2)')  Oml0
      read  ( unit = 30 , fmt = '(70x,   e10.2)')  Omb0
      read  ( unit = 30 , fmt = '(70x,   e10.2)')  hubble
c
      close ( 30 ) 
c
      return
      end

