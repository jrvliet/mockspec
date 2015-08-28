c     --------------
c     Tree structure
c     --------------

      include 'a_numbers.h'
      include 'a_setup.h' 

c     ndim   - number of dimensions
c     nchild - number of children    ( nchild = 2**ndim        )
c     nneigh - number of neighbours  ( nneigh = 2*ndim         )
c     mcell  - number of tree cells  ( mcell  = ncell0 + nchild * moct  )
c     moct   - length of oct arrays  ( moct   = (mcell-ncell0) / nchild )
c     nctot  - total number of cells

c     iOctLv :    >0   - level of an oct
c     iOctPr :         - parent of an oct
c     iOctCh :    >0   - pointer to an oct of children
c                  0   - there are no children; the cell is a leaf
c     iOctNb :    >0   - pointers to neighbouring cells 
c     iOctPs :         - coordinates of Oct centers
c
      parameter ( iFreeLevel = -1000 )    ! fake level for free (unused) cells
c
c     ............positions.................
c
      parameter ( MaxL1     = MaxLevel + 1           )
      parameter ( isize     = 2**MaxL1               )		  
      parameter ( d_x       = 1.0 / 2**MaxL1         ) ! differs from Lesha!
      parameter ( nshift    = nchild - 1             )
      parameter ( nbshift   = nshift - ncell0        ) ! big shift; see Tools
      parameter ( mbshift   = -nbshift               )
      parameter ( ncell01   = -ncell0 - 1            )
c
c     ..........input parameters..............
c
      logical   ladin
      integer nadin
      parameter ( nadin = 100 )

      common /INP01/ radin(nadin), iadin(nadin), ladin(nadin)
c
c     ..........oct information..............
c
      common /TREE01/ iOctPs(ndim,0:moct)       ! coordinates of Oct centers 
      common /TREE02/ iOctNb(nneigh,0:moct)     ! neighbouring cells
      common /TREE03/ iOctPr(0:moct)            ! parents/linked list index
      common /TREE04/ iOctLv(0:moct)            ! Level 
      common /TREE05/ iOctCh(0:nctot)           ! children Octs
      common /TREE06/ iOctFree                  ! first free Oct
      common /TREE07/ nOct                      ! number of Octs in use
      common /TREE08/ iOctMax                   ! the used oct with max. index
c
c     ........linked list of octs............
c
      common /LIST01/ iOctLL1(0:moct)              ! doubly linked list of octs
      common /LIST02/ iOctLL2(0:moct)
      common /LIST03/ iHOLL(MinLevel:MaxLevel+1)  ! linked list header
      common /LIST04/ iNOLL(MinLevel:MaxLevel+1)  ! # of ll entries at a Level
      common /LIST05/ iOLB (2,MinLevel:MaxLevel)  ! boundaries of Level after Sync
c
c     ........physical variables.............
c
c     grav. variables :
c       var (1,*) - total density 
c       var (2,*) - potential (new)
c       var (3,*) - potential (old)
c     hydro variables :
c       hvar(1,*) - gas density 
c       hvar(2,*) - gas energy 
c       hvar(3,*) - x-momentum 
c       hvar(4,*) - y-momentum
c       hvar(5,*) - z-momentum
c       hvar(6,*) - pressure
c       hvar(7,*) - Gamma
c       hvar(8,*) - internal energy 
c
c     during refinement modifications var is used 
c     to store refinement/derefinement indicators
c
      real var(nvar,nctot)
      real hvar(nhvar,nctot)
      real vnw(nhvar-2,nctot)
c
c     ..........common arrays................
c
      common /COMM01/ iCL(nctot)                ! headers of particle LList
      common /VAR01/ var ! vectors of physical variables
      common /VAR02/ hvar
      common /VAR03/ vnw 
      common /VAR04/ gacc(nclmax)
      common /VAR05/ ref(nctot)
      common /VAR07/ RadPre(nctot)
c
c     ..........scratch arrays...............
c
      common /TEMP01/ iSelect(moct) 
      common /TEMP02/ itmp1(moct)
      common /TEMP03/ itmp2(moct)
      common /TEMP04/ itmp3(moct)
      common /TEMP05/ ind(nctot)
      common /TEMP06/ rhophi(2,nclmax)
c
c     .............particle arrays.............
c
c     particles move with different time steps the current 
c     time moment at which particle coordinates are positioned 
c     and its current time step are stored in pt and pdt 
c
      integer np 
      real pt(npmax), pdt(npmax)
      double precision x(npmax), y(npmax), z(npmax), 
     &                vx(npmax), vy(npmax), vz(npmax)
c      real acx(npmax), acy(npmax), acz(npmax)

      common / PART00 /     np 
      common / PART01 /      x
      common / PART02 /      y
      common / PART03 /      z
      common / PART04 /     vx
      common / PART05 /     vy
      common / PART06 /     vz
      common / PART07 / pdummy(npmax)                 ! scratch (for epot)
      common / PART08 / ddummy(npmax)                 ! scratch (for HF)
      common / PART09 /     pt, pdt                   ! particle time & step
      common / PART10 /  iLL(npmax,2)                 ! doubly linked list
      common / PART11 /     acx                       ! particle acceleration
      common / PART12 /     acy 
      common / PART13 /     acz 
      common / PART14 /     pw(npmax)

c
c....  star particles
c
      real tbirth(nstarmax), pw0(nstarmax)
      real zstII(nstarmax), zstIa(nstarmax)
      common / STAR01 / tbirth    ! time of * particles birth in code units
      common / STAR02 / pw0       ! initial mass of stellar particle (current mass is in pw)
      common / STAR03 / zstII, zstIa ! metallicity of the gas the star was born from 

c     .........information about species........
      
      common / SPEC01 /  wpar(nspec)                   ! particle weight
      common / SPEC02 /    sw(nspec)                   ! relative weight
      common / SPEC03 / ekin_(nspec)                   ! kinetic energy 
      common / SPEC04 /  nsp(nspec,2), lsp(nspec)

c     ..........I/O auxiliary arrays..........

      real xpar(npage),ypar(npage),zpar(npage),
     &     vxx(npage),vyy(npage),vzz(npage)

      common / ROW /   xpar,
     &                 ypar,
     &                 zpar,
     &		        vxx, 
     &                  vyy, 
     &                  vzz

      real             recdat(nrecl)                   ! record storage
      equivalence      (recdat(1),xpar(1))              ! position at row

c
c     ..........miscellaneous................
c
      common /MISC01/ CellSize(MinLevel:MaxLevel), 
     &                CellSizei(MinLevel:MaxLevel), 
     &                CellVol (MinLevel:MaxLevel), 
     &                CellVoli(MinLevel:MaxLevel)

c
c     ..... seed for the random number generator ....
c

      integer mrand
      common /RAND01 / mrand
