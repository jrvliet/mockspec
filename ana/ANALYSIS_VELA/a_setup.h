c     ---------------------------------------------
c     code setup : control constants and array sizes
c     ---------------------------------------------

c
c.... parameters for parallelization services
c 
      integer NCPUs , Nproc , LBuffer , MLChunk 
c 
      parameter ( NCPUs    = 24         ) ! anticipated number of processors
      parameter ( Nproc    = NCPUs      )                                 
      parameter ( LBuffer  = 500 000    ) ! Length of buffer for Dens  15000000
      parameter ( MLChunk  = 160 000      ) ! Max Cells per processor for Dens
      parameter ( iSwChunk = 16384      ) ! chunks to work on for Sweep
c 
c.... global setup constants
c 
      integer ndim , mcell , nrow , ng , MinLevel , MaxLevel 
      integer ioff , nchem , nvar , nhvar , nspec
c 
      parameter ( ndim     = 3          ) ! # of spatial dimensions
      parameter ( mcell    = 80 000 000  ) ! # of refinement cells  10000000 (keep <65M cells for CL6csf1)  
      parameter ( nrow     = 128        ) ! # of particles in 1d
      parameter ( nrowreal = 4096      ) ! # real NROW for MM simulations 1024
      parameter ( ng       = 128        ) ! # of 0 lv cells in 1d
      parameter ( MinLevel = 0          ) ! minimum allowed level
      parameter ( MaxLevel = 12          ) ! maximum allowed level
      parameter ( ioff     = 0          ) ! buffer size for arrays
      parameter ( nchem    = 2          ) ! # of chemical species
      parameter ( nMarkers= 0) ! # of non-advected species ( CEVERINO 02032006: The last nMarkers species are non-advected)
      parameter ( nvar     = 3          ) ! # of grav. variables(rho,phi1,phi2)
      parameter ( nhvar    = 8 + nchem  ) ! # of hydro variables
      parameter (nhvarA = nhvar - nMarkers) ! # of advected species . CEVERINO 02032006
      parameter ( nspec    = 6          ) ! # of particle species 1, 2, 3
      parameter ( nsteprun = 1        ) ! Number of steps per each HART run
c 
c     secondary setup constants
c     
c     ..........particles...............
c

CEVERINO03122009      integer npmax , ncold , npage , nrecl
      integer npmax, npage, nrecl

c
  
      parameter ( npmax = 50 000 000 ) 
     

      integer nstarmax 
      parameter ( nstarmax = 15 000 000 )	
c      parameter ( nstarmax =  400 000 )

ceverino03172009      parameter ( ncold     = nrow**3       ) ! # of CDM particles
      parameter ( npage     = nrowreal**2       ) ! # of particles in a page
      parameter ( nrecl     = npage * 6     ) ! length of particle row in words

c
c     .........zero level...............
c
      integer ncell0 , ng2 , narr , nf67 ,  nctot
c
      parameter ( ncell0 = ng**3        ) ! # of zero level cells
      parameter ( nclmax = mcell       ) ! # max # of level cells  (>=ncell0)
      parameter ( ng2    = ng**2        ) ! # of cells in a grid layer
      parameter ( narr   = ng + 1       ) ! FFT parameter
      parameter ( nf67   = ng/2         ) ! FFT parameter
      double precision xn , yn
      parameter ( xn     = 1.d0 + 1.d0*ng) ! boundaries
      parameter ( yn     = 1.d0*ng - 1.0d-6    ) 
      parameter ( nctot  = mcell ) ! total number of cells
c
c     ............tree..................
c
      integer nneigh , neighb , nchild , moct
c
      parameter ( nneigh = 2*ndim       ) ! # of neighbors
      parameter ( neighb = 2*ndim       ) ! # of neighbors
      parameter ( nchild = 2**ndim      ) ! # of children
      parameter ( moct   = (mcell-ncell0)/nchild ) ! # of octs (old nky)

