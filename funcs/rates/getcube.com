c     GLOBAL INFORMATION ABOUT DATA CUBE AND CELL DATA

      integer          ncells            ! number of cells
      integer          Nclevel           ! No. of cells in a level
      integer          Nclevels          ! No. of cell levels
      integer          ID1               ! ID No. of 1st cell in level
      integer          IDN               ! ID No. of last cell in level

      double precision cellsize          ! cell size (pc)
      double precision xposn,yposn,zposn ! cell coordinates (kpc)
      double precision vx,vy,vz          ! cell velocity components (km/s)
      double precision den               ! cell density  (H atoms/cm^3)
      double precision Tcell             ! cell temperature (K)
      double precision SNII,SNIa         ! mass fractions of metals
      double precision clevel

      COMMON/cubedata/ cellsize(maxcells),clevel(maxcells),
     &                 xposn(maxcells),yposn(maxcells),zposn(maxcells),
     &                 vx(maxcells),vy(maxcells),vz(maxcells),
     &                 den(maxcells),Tcell(maxcells),
     &                 SNII(maxcells),SNIa(maxcells)

      COMMON/cubesize/ ncells,Nclevels,Nclevel(Nlmax),
     &                 ID1(Nlmax),IDN(Nlmax)


      double precision Lcell        ! cell size for given cell level
      double precision nHmax,nHmin  ! max,min nH in given cell level
      double precision Tmax,Tmin    ! max,min T in level 
      double precision Zmax,Zmin    ! max,min metal mass frac in level
      double precision Rmax,Rmin    ! max,min galactocentric distance in level

      COMMON/rclevels/ Lcell(Nlmax),
     &                 nHmax(Nlmax),nHmin(Nlmax),
     &                 Tmax(Nlmax),Tmin(Nlmax),
     &                 Zmax(Nlmax),Zmin(Nlmax),
     &                 Rmax(Nlmax),Rmin(Nlmax)
