c     =====================================================================
c                                                                         .
c              ART Version 4:  Fully Threaded Tree Tools                  .
c                                                                         .
c                original algorithms by A.M. Khokhlov                     .                
c                                                                         .
c     =====================================================================
c     -------------------
      function iPr ( iC )      
c     -------------------
c
c     purpose : finds parent for cell iC
c     Input   : iC     - pointers to cells (iC<=ncell0 --> iPr=0)
c
      include 'a_tree.h'
      integer iC
      if ( iC .le. ncell0 ) then 
        iPr = nil 
      else
        iPr = iOctPr(ishft ( iC + nbshift , - ndim ))
      endif
      return
      end

c     -------------------
      function iLv ( iC )      ! 
c     -------------------
c
c     purpose : finds level of cell iC
c     Input   : iC     - pointers to cells (iC<=ncell0 --> iLv=0)

      include 'a_tree.h'
      integer iC
      if ( iC .le. ncell0 ) then 
        iLv = MinLevel
      else
        iLv = iOctLv(ishft ( iC + nbshift , - ndim ))
      endif
      return
      end

c     -------------------------------------
      subroutine Ps ( iC, Posx, Posy, Posz )      ! finds position
c     -------------------------------------
c     purpose : returns coordinates of cell center
c     Input   : iC     - pointer to a cell 
c     Output  : Posxyz    - x,y,z positions
c
      include 'a_tree.h'
      integer iC
      integer idelta(8,3) 
      real Posx, Posy, Posz
      data      idelta / -1,  1, -1,  1, -1,  1, -1,  1,
     &                   -1, -1,  1,  1, -1, -1,  1,  1, 
     &                   -1, -1, -1, -1,  1,  1,  1,  1 /
      if ( iC .le. ncell0 ) then 
        iC_ = iC - 1 
        i = iC_ / ng2               ! ng2 = ng**2 
        j = ( iC_ - i*ng2 ) / ng 
        k =   iC_ - i*ng2 - j*ng
        Posx = float(i) + 1.5
        Posy = float(j) + 1.5
        Posz = float(k) + 1.5
      else
        iC_ = iC + nbshift
        iO = ishft ( iC_ , - ndim )
        id = ishft ( 1, MaxLevel - iOctLv(iO) )   
        j  = iC_ + 1 - ishft( iO , ndim )
        Posx   = d_x * (iOctPs(1,iO) + sign ( id , idelta(j,1) ))
        Posy   = d_x * (iOctPs(2,iO) + sign ( id , idelta(j,2) ))
        Posz   = d_x * (iOctPs(3,iO) + sign ( id , idelta(j,3) )) 
      endif
      return
      end

c     -----------------------
      function iCh ( iC , j )
c     -----------------------
c
c     purpose : finds j-th child
c     Input   : iC     - pointers to cells (iC<=0 --> iCh=0)
c               j      - child ( 1 - 8 )
c
      include 'a_tree.h'
      integer iC
      iOC = iOctCh(iC)
      if ( iOC .le. nil ) then
        iCh = nil
      else
        iCh = ( iOC - 1 ) * nchild + j + ncell0
      end if
      return
      end

c     ---------------------------------
      subroutine iChAll ( iC , iChild )
c     ---------------------------------
c
c     purpose : finds all children, iChild(nchild) of cell iC
c     Input   : iC     - pointer of the parent cells (iC<=0 --> iCh=0)
c     Output  : iChild - array of nchild pointers to the children cells
c
      include 'a_tree.h'
      integer iC, iChild(nchild), j, i0
c
      iOC = iOctCh(iC)
      if ( iOC .le. nil ) then
        iCh = nil
        iChild(1) = nil
        return
      else
        i0 = ( iOC - 1 ) * nchild + ncell0
        do j = 1 , nchild
          iChild(j) = i0 + j
        enddo
      end if
      return
      end

c     -------------------------------
      function iNb ( iC , idir )      
c     -------------------------------

c     Purpose : finds neighbour of cell iC in idir direction
c     Input   : iC     - pointer to cells 
c               idir   - direction ( 1 - 6 )
c
      include 'a_tree.h'
      integer iC, idir
      parameter ( log2ng = 7 , log2ng2 = log2ng*2 , ng1 = ng - 1) 
      dimension lnbr(0:8,6), kiot(0:8,6), ishift(6,3)
      data lnbr / 0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 5, 6, 7, 8, 1, 2, 3, 4,      ! idir = 5,6
     &            0, 5, 6, 7, 8, 1, 2, 3, 4/      ! idir = 5,6
      data kiot / 1, 0, 1, 0, 1, 0, 1, 0, 1,      ! idir = 1
     &            1, 1, 0, 1, 0, 1, 0, 1, 0,      ! idir = 2
     &            1, 0, 0, 1, 1, 0, 0, 1, 1,      ! idir = 3
     &            1, 1, 1, 0, 0, 1, 1, 0, 0,      ! idir = 4
     &            1, 0, 0, 0, 0, 1, 1, 1, 1,      ! idir = 5
     &            1, 1, 1, 1, 1, 0, 0, 0, 0/      ! idir = 6
      data ishift / -1, 1,  0, 0,  0, 0,
     &               0, 0, -1, 1,  0, 0,
     &               0, 0,  0, 0, -1, 1  /
c      Return
      if ( iC .le. ncell0 ) then 
        i = ( iC - 1 ) / ng2            
        j = ( iC - 1 - i*ng2  ) / ng 
        k =   iC - i*ng2 - j*ng    
        i = i + ishift(idir,1)
        j = j + ishift(idir,2)
        k = k + ishift(idir,3)
        if ( i .gt. ng1 ) i = 0              ! periodic box check
        if ( j .gt. ng1 ) j = 0
        if ( k .gt. ng ) k = 1
        if ( i .lt.  0 ) i = ng1
        if ( j .lt.  0 ) j = ng1
        if ( k .lt.  1 ) k = ng
        iNb = (i*ng + j)*ng + k
c        iNb = ishft(i,log2ng2) + ishft(j,log2ng) + k
      else

        iO    = ishft ( iC + nbshift , -ndim )  
        iO_   = iO -1                           ! cells oct
        icloc = iC -ncell0 -ishft( iO_ , ndim ) ! cells local # ( 1 - 8 )
        inl   = lnbr(icloc,idir) +ncell0        ! neighbours local number
        If(kiot(icloc,idir).eq.0)Then          ! 0/1 - out/in key
           iON   = iOctNb(idir,iO)             ! cells octs neighbouring cell
           iONC  = iOctCh(iON)                 ! neighbouring oct 
           If(iONC.eq.0)Then
              iNb   =  iON
           Else
              iNb   =  inl  +ishft(iONC-1,ndim)
           EndIf
        Else
              iNb   =  inl  +ishft ( iO_, ndim ) 
        EndIf
      end if
      return
      end


c     -------------------------------
      subroutine NbAll0 ( iC , iNbb )
c     -------------------------------
      include 'a_tree.h'
      parameter ( ngg = ng2 + ng )
      integer iC, iNbb(6)

      i = ( iC - 1 ) / ng2 + 1
      j = ( iC - 1 - (i-1)*ng2  ) / ng + 1
      kngg =  iC - i*ng2 - j*ng 
      k =   kngg + ngg
      i1 = i -1                           ! 1st direction
      if ( i1 .lt.  1 ) i1 = ng
      i2 = i +1                           ! 2nd direction
      if ( i2 .gt. ng ) i2 = 1
      j1 = j -1                           ! 3rd direction
      if ( j1 .lt.  1 ) j1 = ng
      j2 = j +1                           ! 4th direction
      if ( j2 .gt. ng ) j2 = 1
      k1 = k -1                           ! 5rd direction
      if ( k1 .lt.  1 ) k1 = ng
      k2 = k +1                           ! 6th direction
      if ( k2 .gt. ng ) k2 = 1
      
      jng  = j*ng
      jngg = jng - ngg
      ing2 = i*ng2

      iNbb(1) = i1*ng2 + jng  + kngg
      iNbb(2) = i2*ng2 + jng  + kngg
      iNbb(3) = j1*ng  + ing2 + kngg
      iNbb(4) = j2*ng  + ing2 + kngg
      iNbb(5) =  k1    + ing2 + jngg
      iNbb(6) =  k2    + ing2 + jngg

      return
      end

c     -------------------------------
      subroutine NbAllL ( iC , iNbb )
c     -------------------------------
c
c     select all 6 neighbors of cell iC of level > 0
c
      include 'a_tree.h'
      integer iC, iNbb(6)
      integer iINb(8,3)   ! internal sibling neighbors number
      integer iEDir(8,3)  ! directions to external neighbors
      integer iIDir(8,3)  ! directions to internal neighbors
      integer lnbr(8,6)

      data iINb  / 2, 1, 4, 3, 6, 5, 8, 7,     
     &             3, 4, 1, 2, 7, 8, 5, 6,
     &             5, 6, 7, 8, 1, 2, 3, 4  / 

      data iEDir / 1, 2, 1, 2, 1, 2, 1, 2,
     &             3, 3, 4, 4, 3, 3, 4, 4,
     &             5, 5, 5, 5, 6, 6, 6, 6  /

      data iIDir / 2, 1, 2, 1, 2, 1, 2, 1,
     &             4, 4, 3, 3, 4, 4, 3, 3,
     &             6, 6, 6, 6, 5, 5, 5, 5  /

      data lnbr / 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            5, 6, 7, 8, 1, 2, 3, 4,      ! idir = 5,6
     &            5, 6, 7, 8, 1, 2, 3, 4/      ! idir = 5,6

      iO  = ishft ( iC + nbshift , -ndim )
      iO_ = iO - 1                              ! cells oct
      ibase = ncell0 + ishft( iO_ , ndim )
      icloc = iC - ibase                        ! cells local # ( 1 - 8 )
      
      iNbb(iIDir(icloc,1)) = iINb(icloc,1) + ibase 
      iNbb(iIDir(icloc,2)) = iINb(icloc,2) + ibase 
      iNbb(iIDir(icloc,3)) = iINb(icloc,3) + ibase 

      idir  = iEDir(icloc,1)
      iON   = iOctNb(idir,iO)             ! cells octs neighbouring cell
      iONC  = iOctCh(iON)                 ! neighbouring oct
      if ( iONC .eq. 0 ) then
        iNbb(idir) =  iON
      else
        iNbb(idir) =  lnbr(icloc,idir) + ncell0 + ishft(iONC-1,ndim)
      endif

      idir  = iEDir(icloc,2)
      iON   = iOctNb(idir,iO)             ! cells octs neighbouring cell
      iONC  = iOctCh(iON)                 ! neighbouring oct
      if ( iONC .eq. 0 ) then
        iNbb(idir) =  iON
      else
        iNbb(idir) =  lnbr(icloc,idir) + ncell0 + ishft(iONC-1,ndim)
      endif

      idir  = iEDir(icloc,3)
      iON   = iOctNb(idir,iO)             ! cells octs neighbouring cell
      iONC  = iOctCh(iON)                 ! neighbouring oct
      if ( iONC .eq. 0 ) then
        iNbb(idir) =  iON
      else
        iNbb(idir) =  lnbr(icloc,idir) + ncell0 + ishft(iONC-1,ndim)
      endIf
      
      return
      end

c     -------------------------------
      subroutine iNbAll ( iC , iNbb )      
c     -------------------------------

c     Purpose : finds neighbour of cell iC in idir direction
c     Input   : iC     - pointer to cells 
c               idir   - direction ( 1 - 6 )
c
      include 'a_tree.h'
      parameter ( ngg = ng2 + ng )
      integer iC, idir
      dimension lnbr(0:8,6), kiot(0:8,6), ishift(6,3),
     &          iNbb(6) 
      data lnbr / 0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 2, 1, 4, 3, 6, 5, 8, 7,      ! idir = 1,2
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 3, 4, 1, 2, 7, 8, 5, 6,      ! idir = 3,4
     &            0, 5, 6, 7, 8, 1, 2, 3, 4,      ! idir = 5,6
     &            0, 5, 6, 7, 8, 1, 2, 3, 4/      ! idir = 5,6
      data kiot / 1, 0, 1, 0, 1, 0, 1, 0, 1,      ! idir = 1
     &            1, 1, 0, 1, 0, 1, 0, 1, 0,      ! idir = 2
     &            1, 0, 0, 1, 1, 0, 0, 1, 1,      ! idir = 3
     &            1, 1, 1, 0, 0, 1, 1, 0, 0,      ! idir = 4
     &            1, 0, 0, 0, 0, 1, 1, 1, 1,      ! idir = 5
     &            1, 1, 1, 1, 1, 0, 0, 0, 0/      ! idir = 6
      data ishift / -1, 1,  0, 0,  0, 0,
     &               0, 0, -1, 1,  0, 0,
     &               0, 0,  0, 0, -1, 1  /
c      Return
      if ( iC .le. ncell0 ) then 
        i = ( iC - 1 ) / ng2 + 1            
        j = ( iC - 1 - (i-1)*ng2  ) / ng + 1
        k =   iC - i*ng2 - j*ng + ngg    
        i1 = i -1                           ! 1st direction
        if ( i1 .lt.  1 ) i1 = ng
        i2 = i +1                           ! 2nd direction
        if ( i2 .gt. ng ) i2 = 1            
        j1 = j -1                           ! 3rd direction
        if ( j1 .lt.  1 ) j1 = ng
        j2 = j +1                           ! 4th direction
        if ( j2 .gt. ng ) j2 = 1            
        k1 = k -1                           ! 5rd direction
        if ( k1 .lt.  1 ) k1 = ng
        k2 = k +1                           ! 6th direction
        if ( k2 .gt. ng ) k2 = 1            

        iNbb(1) = i1*ng2 + j *ng  + k  -ngg
        iNbb(2) = i2*ng2 + j *ng  + k  -ngg
        iNbb(3) = j1*ng  + i *ng2 + k  -ngg
        iNbb(4) = j2*ng  + i *ng2 + k  -ngg
        iNbb(5) =  k1    + i *ng2 + j *ng  -ngg
        iNbb(6) =  k2    + i *ng2 + j *ng  -ngg       
      else

        iO    = ishft ( iC + nbshift , -ndim )  
        iO_   = iO -1                           ! cells oct
        icloc = iC -ncell0 -ishft( iO_ , ndim ) ! cells local # ( 1 - 8 )
        Do idir =1,6
          If(kiot(icloc,idir).eq.0)Then          ! 0/1 - out/in key
           iON   = iOctNb(idir,iO)             ! cells octs neighbouring cell
           iONC  = iOctCh(iON)                 ! neighbouring oct 
           If(iONC.eq.0)Then
              iNbb(idir) =  iON
           Else
              inl        =  lnbr(icloc,idir) +ncell0 ! neighbours local number
              iNbb(idir) =  inl +ishft(iONC-1,ndim)
           EndIf
          Else
              inl        =  lnbr(icloc,idir) +ncell0 ! neighbours local number
              iNbb(idir) =  inl +ishft ( iO_, ndim ) 
          EndIf
        EndDo
      end if
      return
      end

c     ----------------------------------------
      function iFindCell0 ( L , xp , yp , zp )
c     ----------------------------------------
c
c     purpose: finds a cell of Level enclosing 
c              point with coordinates xp, yp, zp
c    
c     input  : Level, xp, yp, zp
c
c     this routine differs from iFindCell in that it returns 0
c     if only cell of level < L can be found around xp, yp, zp
c
      include 'a_tree.h'
      parameter ( ng1 = ng + 1 )
      integer iCell, L
      double precision xp, yp, zp, xx, yy, zz

      ip    = int(xp)
      jp    = int(yp)
      kp    = int(zp)
      if ( ip .eq. 0 ) ip = ng 
      if ( jp .eq. 0 ) jp = ng
      if ( kp .eq. 0 ) kp = ng
      if ( ip .eq. ng1 ) ip = 1
      if ( jp .eq. ng1 ) jp = 1
      if ( kp .eq. ng1 ) kp = 1

      iCell = iCellIndex ( ip , jp , kp )

      if ( L .eq. MinLevel ) then
        iFindCell0 = iCell
        return
      endif

      xx    = xp - ip
      yy    = yp - jp
      zz    = zp - kp

      iter  = 0
 1    continue
      dxr   = ishft(1,iter+1)  ! = 2**iter
      iChild = 1 +     mod ( int ( xx * dxr ) , 2 )
     &           + 2 * mod ( int ( yy * dxr ) , 2 )
     &           + 4 * mod ( int ( zz * dxr ) , 2 )
      iCCh = iOctCh(iCell)
      if ( iCCh .gt. nil ) then
        iter = iter + 1
        iCell = ncell0 + ishft ( iCCh - 1 , ndim ) + iChild 
      endif
      if ( (iter.lt.L) .and. (iOctCh(iCell).gt.nil) ) go to 1

      if ( iter .eq. L ) then 
        iFindCell0 = iCell
      else
        iFindCell0 = nil 
      endif           

      return
      end

c     ---------------------------------------
      function iFindCell ( L , xp , yp , zp )
c     ---------------------------------------
c
c     purpose: finds a cell of Level or above enclosing 
c              point with coordinates xp, yp, zp
c     input  : Level, xp, yp, zp
c
      include 'a_tree.h'
      parameter ( ng1 = ng + 1 )
      integer iCell, L
      double precision xp, yp, zp, xx, yy, zz

      ip    = int(xp)
      jp    = int(yp)
      kp    = int(zp)
      if ( ip .eq. 0 ) ip = ng 
      if ( jp .eq. 0 ) jp = ng
      if ( kp .eq. 0 ) kp = ng
      if ( ip .eq. ng1 ) ip = 1
      if ( jp .eq. ng1 ) jp = 1
      if ( kp .eq. ng1 ) kp = 1

      iCell = iCellIndex ( ip , jp , kp )

      if ( L .eq. MinLevel ) then
        iFindCell = iCell
        return
      endif

      xx    = xp - ip
      yy    = yp - jp
      zz    = zp - kp

      iter  = 0
 1    iter  = iter + 1
      dxr   = ishft(1,iter)  ! = 2**iter
      iChild = 1 +     mod ( int ( xx * dxr ) , 2 )
     &           + 2 * mod ( int ( yy * dxr ) , 2 )
     &           + 4 * mod ( int ( zz * dxr ) , 2 )
      iCCh = iOctCh(iCell)
      if ( iCCh .gt. nil ) then 
        iCell = ncell0 + ishft ( iCCh - 1 , ndim ) + iChild 
      endif
      if ( (iter.lt.L) .and. (iOctCh(iCell).gt.nil) ) go to 1

      iFindCell = iCell
            
      return
      end

c     -------------------
      function iGetOct ()       ! returns free Oct, or zero
c     -------------------
      include 'a_tree.h'
!MIC$GUARD
      i = iOctFree
      if ( iOctFree .gt. 0 ) then
        iOctFree = iOctPr(i)
        nOct = nOct + 1
        iOctPr(i) = 0
        do k = 1 , ndim
          iOctPs(k,i) = 0
        end do
        do k = 1 , nneigh
          iOctNb(k,i) = 0
        end do
        iOctLL1(i) = 0
        iOctLl2(i) = 0
        do k = 1 , nchild
          iOctCh(ncell0+(i-1)*nchild+k) = 0
        end do
      else
        if ( iOctFree .lt. 0 ) then 
          write(*,*) 'iGetOct error: iOctFree =',iOctFree
          stop
        endif
      end if
      iGetOct = i
!MIC$END GUARD
      return
      end

c     -----------------------
      function iFreeOct ( i )   ! returns i to the pool of free Octs
c     -----------------------
      include 'a_tree.h'
!MIC$GUARD
      iOctPr(i)  = iOctFree
      iOctLv(i)  = iFreeLevel          ! oct = free
      iOctLL1(i) = 0
      iOctLL2(i) = 0
      iOctFree  = i
      nOct      = nOct - 1
      iFreeOct  = 0
!MIC$END GUARD
      return
      end

c     ------------------------------------
      function iCellIndex ( ic , jc , kc )  
c     ------------------------------------
c
c     returns index of parent cell from iPr 
c     which corresponds to given combination of i,j,k
c
c     (!) Note: algorithm of cell identification depends on the
c     order in which we assume cells are placed in iPr:
c     assumed order is  i,j,k. When k is changing most rapidly
c     (see make_root_mesh subroutine)
c
      include 'a_setup.h'

      iCellIndex = ((ic-1)*ng + (jc-1))*ng + kc

      return
      end

c     --------------------------------------------
      subroutine iCellijk ( icell , ic , jc , kc )
c     --------------------------------------------
c
c     for a zero-level cell returns 
c     corresponding i,j,k indices
c
      include 'a_setup.h'

      ic = (icell - 1)/ng2 + 1
      jc = (icell - (ic-1)*ng2 - 1)/ng + 1
      kc =  icell - (ic-1)*ng2 - (jc-1)*ng

      return
      end

c     -----------------------------
      subroutine Get_MaxLevelNow ()
c     -----------------------------
c
      include 'a_tree.h'
      include 'a_control.h'

      MaxLevelNow = MinLevel
      if ( iNOLL(MinLevel+1) .eq. nil ) then 
        MaxLevelNow = MinLevel
        return
      else
        do Level = MaxLevel , MinLevel + 1 , -1
          if ( iNOLL(Level) .gt. nil ) then 
            MaxLevelNow = Level 
            return
          endif
        enddo
      endif
      MaxLevelNow = Level - 1
      return
      end

c     ------------------------------
      subroutine iKy_Insert ( iOct )
c     ------------------------------
c
c     purpose: inserts iOct in doubly linked list iOctLL1 , iOctLL2
c    
      include 'a_tree.h'
      integer iOct

      Level = iOctLv(iOct)
      iOctLL1(iOct) = iHOLL(Level)

      if ( iHOLL(Level) .ne. nil ) then    ! if Level was not empty
        iOctLL2(iHOLL(Level)) = iOct
      endif
      iHOLL(Level) = iOct
      iNOLL(Level) = iNOLL(Level) + 1
      iOctLL2(iOct) = nil

      return
      end

c     ------------------------------
      subroutine iKy_Delete ( iOct )
c     ------------------------------
c
c     deletes iOct from the linked list of octs
c
      include 'a_tree.h'
      integer iOct

      Level = iOctLv(iOct)

      if ( iOctLL2(iOct) .ne. nil ) then
        iOctLL1(iOctLL2(iOct)) = iOctLL1(iOct)
      else
        iHOLL(Level) = iOctLL1(iOct)
      endif

      if ( iOctLL1(iOct) .ne. nil ) then
        iOctLL2(iOctLL1(iOct)) = iOctLL2(iOct)
      endif

      iNOLL(Level) = iNOLL(Level) - 1
      if ( iNOLL(Level) .eq. nil ) then 
        iHOLL(Level) = nil
      endif

      return
      end

c     ------------------------------------------
      integer function iGet_Oct_Number ( icell ) 
c     ------------------------------------------
c
c     purpose: auxiliary function to map real index 
c              icell onto (icell-ncell0)/nchild space
c     note   : icell should be of Level > MinLevel 
c
      include 'a_tree.h'
      integer icell
      iGet_Oct_Number = ishft ( icell + nbshift , - ndim )
      return
      end

c     ---------------------------------
      function iGet_Cell_Index ( iOct ) 
c     ---------------------------------
c
c     returns index of the first cell in iOct
c
      include 'a_tree.h'
      integer iOct
      iGet_Cell_Index = (iOct - 1) * nchild + ncell0 + 1
      return
      end

c     -----------------------------------------
      subroutine Select_Cells ( Level , nLevel ) 
c     -----------------------------------------
c
c     purpose: selects cells (every eights) of a specified Level
c
      include 'a_tree.h'
      include 'a_control.h'

      integer Level , nLevel 

      nLevel = iNOLL(Level)     ! get level boundaries for index array 

      icell  = iGet_Cell_Index ( iHOLL(Level) )
      do j = 1 , nLevel         
        iSelect(j) = icell
        iOct       = (icell - ncell0) / nchild + 1
        iCell      = (iOctLL1(iOct) - 1) * nchild + ncell0 + 1
      enddo
      return
      end

c     ---------------------------
      function iChildNumber ( i )
c     ---------------------------
c
c     returns the child number 
c     for a cell i; should be iLv(i) > MinLevel
c
      include 'a_tree.h'
      if ( iLv(i) .le. MinLevel ) then
        write(*,*) 'Error in iChildNumber: MinLevel cell'
        stop
      endif
      iChildNumber = mod ( i - ncell0 - 1 , nchild ) + 1
      return
      end

c     -------------------------------------------------
      function iChildNumber2 ( xint , yint , zint , i )
c     -------------------------------------------------
c
c     returns a childnumber i a parent i 
c     to which point x,y,z is belonging
c
      include 'a_tree.h'
      real Posx, Posy, Posz

      call Ps ( i , Posx, Posy, Posz )

      ifl = sign ( one , (xint - Posx) )
      jfl = sign ( one , (yint - Posy) )
      kfl = sign ( one , (zint - Posz) )
      iChildNumber2 = (ifl + 1)/2 + (jfl + 1) + 2 * (kfl + 1) + 1

      return
      end

c     ------------------------------------
      function iChildNumber3 ( i , j , k )
c     ------------------------------------
c     computes Child # for given ( i , j , k ) triple:
c
c     1 = ( -1 , -1 , -1 )   5 = ( -1 , -1 ,  1 )
c     2 = (  1 , -1 , -1 )   6 = (  1 , -1 ,  1 )
c     3 = ( -1 ,  1 , -1 )   7 = ( -1 ,  1 ,  1 )
c     4 = (  1 ,  1 , -1 )   8 = (  1 ,  1 ,  1 )
c
c     note: all operations ARE integer
c     
      integer i , j , k
      iChildNumber3 = (i + 1)/2 + (j + 1) + 2 * k + 3
      return
      end
