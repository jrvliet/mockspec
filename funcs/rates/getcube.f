C     SUBROUTINE getcube
!     Driver for reading in gas cell information stored in a *GZ* file 
!     and writing it back out.
!
!     Example of a GZ file: MW9_GZ889.a1.001.txt
!
!
!     Author: Liz Klimek
!             Chris Churchill (added cell level data gathering)
!
c     *** IMPORTANT *** the cell data MUST be sorted in order of
c     decreasing cell size of cell level algorithm will fail
c
c     - SORT THE GRID BEFORE RUNNING THE PROGRAM
c
c

!-----------------------------------------------------------------------

      SUBROUTINE readGZ(GZfile)
     
!     Read in a GZ data cube, which contains information about the gas 
!     cells in a simulation.  Store each column of data into an array.
!
!     Example of a GZ file: MW9_GZ889.a1.001.txt

      implicit none
      
      include          'getcube.h'
      integer          i,knt
      double precision R
      character*80     GZfile      
      character*80     headerline

      include          'getcube.com'


      CALL notify(' ',-1)
      CALL notify(' HYDRO DATA CUBE',-1)
      CALL notify(' processing gas cell data...',-1)
     
      
      open(unit=12,file=GZfile)
      
      do i=1,nhdr  ! NHDR defined in getcube.h
        read(12,*) headerline
      enddo

      i            = 1
      clevel(i)    = 1
      knt          = 1
      Nclevel(knt) = 1

      read(12,*) cellsize(i),xposn(i),yposn(i),zposn(i),
     +           vx(i),vy(i),vz(i),den(i),Tcell(i),
     +           SNII(i),SNIa(i)

      Lcell(knt) = cellsize(i)
      nHmax(knt) = den(i)
      nHmin(knt) = den(i)
      Tmax(knt)  = Tcell(i)
      Tmin(knt)  = Tcell(i)
      Zmax(knt)  = SNII(i)+SNIa(i)        
      Zmin(knt)  = Zmax(knt)
      Rmax(knt)  = sqrt(xposn(i)**2+yposn(i)**2+zposn(i)**2)
      Rmin(knt)  = Rmax(knt)

      i = 2
      DO
        READ(12,*,END=6) cellsize(i),xposn(i),yposn(i),zposn(i),
     +                   vx(i),vy(i),vz(i),den(i),Tcell(i),
     +                   SNII(i),SNIa(i)
        IF (cellsize(i).ne.cellsize(i-1)) then  ! increment to next level
          knt          = knt + 1
          clevel(i)    = knt
          Nclevel(knt) = 1
          Lcell(knt)   = cellsize(i)
          nHmax(knt)   = den(i)
          nHmin(knt)   = den(i)
          Tmax(knt)    = Tcell(i)
          Tmin(knt)    = Tcell(i)
          Zmax(knt)    = SNII(i)+SNIa(i)        
          Zmin(knt)    = Zmax(knt)
          Rmax(knt)    = sqrt(xposn(i)**2+yposn(i)**2+zposn(i)**2)
          Rmin(knt)    = Rmax(knt)
        ELSE                                   ! statistics for this level
          clevel(i)    = knt
          Nclevel(knt) = Nclevel(knt) + 1
          nHmax(knt)   = max(nHmax(knt),den(i))
          nHmin(knt)   = min(nHmin(knt),den(i))
          Tmax(knt)    = max(Tmax(knt),Tcell(i))
          Tmin(knt)    = min(Tmin(knt),Tcell(i))
          Zmax(knt)    = max(Zmax(knt),SNII(i)+SNIa(i))
          Zmin(knt)    = min(Zmin(knt),SNII(i)+SNIa(i))
          R            = sqrt(xposn(i)**2+yposn(i)**2+zposn(i)**2)
          Rmax(knt)    = max(Rmax(knt),R)
          Rmin(knt)    = min(Rmin(knt),R)
       END IF  
        i=i+1
        IF (i.gt.maxcells) GOTO 1000
      ENDDO
      
c     close the file

6     CLOSE(12)
      

c     store the total number of cells read in and communicate

      ncells = i-1   

c     error trap 1

      IF (ncells.gt.maxcells) then
       WRITE(6,*) '...................................................'
       WRITE(6,*) 'WARNING: number of cells exceeds physical dimension'
       WRITE(6,*) 'defined in the getcube.h file.  Number of cells in'
       WRITE(6,*) 'current cube is NCELLS= ',ncells
       WRITE(6,*) 'physical dimension is MAXCELLS= ',maxcells
       WRITE(6,*) 'this error results in cell data curruption'
       WRITE(6,*) 'TO FIX- edit the file getcube.h and increase the'
       WRITE(6,*) 'value of parameter MAXCELLS and re-make rates.'
       WRITE(6,*) '...................................................'
       STOP
      END IF

      CALL notify('  No. of grid cells  = ',ncells)

      Nclevels = knt

      CALL notify('  No. of cell levels = ',Nclevels)
      CALL notify('  Cell gas characteristics:',-1)
      CALL notify('  - Lcell, Rmin, Rmax in parsecs',-1)
      CALL notify('  - nHmin, nHmax, Tmin, Tmax in cgs, Kelvin',-1)
      CALL notify('  - Zmin, Zmax absolute mass metal fraction',-1)
      CALL notify('  ----------------------------------------------',-1)

      WRITE(6,599) '  Level','Ncells','ID1','IDN','Lcell','Rmin','Rmax',
     &                'nHmin','nHmax','Tmin','Tmax','Zmin','Zmax'
      WRITE(4,599) '  Level','Ncells','ID1','IDN','Lcell','Rmin','Rmax',
     &                'nHmin','nHmax','Tmin','Tmax','Zmin','Zmax'
      DO knt=1,Nclevels
       IF (knt.eq.1) then
         ID1(knt) = 1 
         IDN(knt) = Nclevel(knt)
       ELSE
         ID1(knt) = ID1(knt-1) + Nclevel(knt-1) 
         IDN(knt) = IDN(knt-1) + Nclevel(knt)
       END IF
       WRITE(6,600) knt,Nclevel(knt),ID1(knt),IDN(knt),Lcell(knt),
     &            Rmin(knt),Rmax(knt),nHmin(knt),nHmax(knt),
     &            Tmin(knt),Tmax(knt),Zmin(knt),Zmax(knt)
       WRITE(4,600) knt,Nclevel(knt),ID1(knt),IDN(knt),Lcell(knt),
     &            Rmin(knt),Rmax(knt),nHmin(knt),nHmax(knt),
     &            Tmin(knt),Tmax(knt),Zmin(knt),Zmax(knt)
      ENDDO
      CALL notify('  ----------------------------------------------',-1)

      RETURN

c     error trap 2

 1000 ncells = i - 1

      WRITE(6,*) '...................................................'
      WRITE(6,*) 'WARNING: number of cells exceeds physical dimension'
      WRITE(6,*) 'defined in the getcube.h file.  Number of cells in'
      WRITE(6,*) 'current cube is NCELLS= ',ncells
      WRITE(6,*) 'physical dimension is MAXCELLS= ',maxcells
      WRITE(6,*) 'this error results in cell data curruption'
      WRITE(6,*) 'TO FIX- edit the file getcube.h and increase the'
      WRITE(6,*) 'value of parameter MAXCELLS and re-make rates.'
      WRITE(6,*) '...................................................'
      STOP

 599  FORMAT(1x,a5,2x,3a8,3x,3a7,10a11)
 600  FORMAT(1x,i5,2x,3i8,3x,3f7.1,1p10e11.3)

      END
      
!-----------------------------------------------------------------------

      SUBROUTINE writeGZcell(i,m,kdo,jdo)
     
!     Write information about gas cells to file.
c     most mods from Chris Churchill are in this subroutine


c     i      = cell index
c     j      = ion file index
c     kdo(M) = atomic number of ion
c     jdo(M) = ion stage of ion
 

      include          'rates.h'
      include          'getcube.h'
      integer          i,m,k,j
      integer          kdo(Imax),jdo(Imax)
      double precision ionden,specden,Lcellkpc,ionfrac,alphsol,abund
      character*80     nlabel,slabel,flabel
      include          'rates.com'
      include          'getcube.com'

c     convert cell size to kiloparsecs on output

      Lcellkpc = 1.0d-3*cellsize(i)  

c     set the ion information and write

       k = kdo(m)
       j = jdo(m)

      IF (k.eq.0) then  ! electron density

       nlabel  = 'ne'
       slabel  = 'ne'
       flabel  = 'fne'
       ionden  = eden
       specden = eden
       ionfrac = 1.0d0
       alphsol = 0.0d0
       abund   = 0.0d0
      
      ELSE  ! ion k,j

       nlabel  = 'n'//ionID(k,j)
       slabel  = 'n'//specID(k)
       flabel  = 'f'//ionID(k,j)
       ionden  = nkj(k,j)
       specden = nk(k)
       ionfrac = fion(k,j)
       alphsol = alphksol(k)
       abund   = fabund(k)

      END IF

      IF (i.eq.1) then 
       write(m+10,200) '1','2','3','4','5','6','7','8',
     +                 '9','10','11','12','13','14','15',
     +                 '16','17','18','19','20','21'
       write(m+10,201) 'cell size','x','y','z',
     +                 'vx','vy','vz',
     +                 'nH','temperature',
     +                 'SNII mass frac','SNIa mass frac',
     +                  slabel,flabel,nlabel,
     +                  'alpha_sol','alpha_Zmet',
     +                  'ID','t_ph','t_rec','t_coll','t_cool'
      ENDIF

      write(m+10,'(es11.4, 15es16.4e3, 4x, i8, 2x, 4es11.3)')
     +        Lcellkpc,xposn(i),yposn(i),zposn(i),
     +        vx(i),vy(i),vz(i),nH,Tcell(i),
     +        SNII(i),SNIa(i),specden,ionfrac,ionden,
     +        alphsol,abund,i,tau_ph(m),tau_rec(m),
     +        tau_coll(m),tau_cool
     
      RETURN

 200  FORMAT(5x,a,9(15x,a),6(14x,a),14x,a,6x,a,4(10x,a))
 201  FORMAT(a10,11x,a1,15x,a1,15x,a1,15x,a2,14x,a2,14x,a2,14x,a2,
     +          9x,a11,3x,a14,3x,a14,7x,a7,9x,a7,9x,a7,
     +          7x,a9,6x,a10,4x,a8,3x,a8,3x,a8,3x,a8,3x,a8)

      END
      
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      SUBROUTINE writecooldat(i,z,T,Zmet)
     
!     Write out box with cooling function data (for post-post prociessing)

c     i      = cell index

      implicit none
      include          'rates.h'
      integer          kprnt
      parameter        (kprnt=11)
      integer          i,ii,len,k(kprnt),kk
      double precision z,T,Zmet
      character*80     slabel(kprnt)
      include          'rates.com'

c     the k indices of the species to print

      data (k(kk),kk=1,11)/1,2,6,7,8,10,12,14,16,20,26/

c     make the header labels for the densities

      DO kk=1,kprnt
       len = 0
       DO ii=1,80
        IF(specID(k(kk))(ii:ii).eq.' ') then
         len = ii-1
         GOTO 01
        ENDIF
       ENDDO
 01    slabel(kk) = 'n'//specID(k(kk))(1:len)//'/nH'
      ENDDO

      IF (i.eq.1) then 
       write(77,200) '1','2','3','4','5','6','7','8',
     +                 '9','10','11','12','13','14','15',
     +                 '16','17','18'
       write(77,201) 'ID','z','nH','T','Zcell','ne','ntot',
     +                (slabel(kk),kk=1,kprnt)
      ENDIF

      write(77,202) i,z,nH,T,Zmet,eden,nions,
     +              ((nk(k(kk))/nH),kk=1,kprnt)

      RETURN

 200  FORMAT(8x,a,6x,a,7(10x,a),9(9x,a))
 201  FORMAT(7x,a,6x,a1,9x,a2,10x,a1,8x,a5,7x,a2,7x,a4,4x,11(5x,a6))
 202  FORMAT(1x,i8,17(1pe11.2))

      END
      
!-----------------------------------------------------------------------
