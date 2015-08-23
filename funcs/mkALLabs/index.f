c     PACKAGE: Mockspec
c     PROGRAM: general module for multiple programs in the package
c     MODULE:  index.f
c
c     DESCRIPTION
c     a modified index sorting routine from Numerical Recipes
c
c
c.............................................................................
c

      SUBROUTINE index(n,arrin,indx)

c     indexes the character array arrin of length n. outputs the array
c     indx such that arrin(index(j)) (j=1,...n) is in ascending
c     order. the input array arrin is not changed.

c     this is a modified Numerical Recipe routine originally called
c     indexx; the modification is that array arrin is changed to a
c     character array from a real array

c     for the cullabs and cellabs programsq; this routine obtains the
c     ascending sorted index order of the ion names (i.e., HI, CIV,
c     MgII) so that we can easily determine the number of unique ion
c     names in the LOS data output from program Mockspec/anal/sysanal

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none

      integer          n,indx(n),indxt,j,l,ir,i
      character*80     arrin(n),q



c     iniitalize the index array with consecutive integers

      do 11 j=1,n
        indx(j) = j
11    continue

c     from here on out its heapsort, but with indirect indexing through
c     indx in all references to arrin. 

      l  = n/2+1
      ir = n

10    continue

        if (l.gt.1) then
          l     = l-1
          indxt = indx(l)
          q     = arrin(indxt)
        else
          indxt    = indx(ir)
          q        = arrin(indxt)
          indx(ir) = indx(1)
          ir       = ir-1
          if (ir.eq.1) then
            indx(1) = indxt
            RETURN         ! the return condition
          endif
        endif

        i = l
        j = l+l

20      if (j.le.ir) then
         if (j.lt.ir) then
            if (arrin(indx(j)).lt.arrin(indx(j+1))) j = j+1
          endif
          if (q.lt.arrin(indx(j))) then
            indx(i) = indx(j)
            i       = j
            j       = j+j
          else
            j = ir+1
          endif
        go to 20
        endif
        indx(i) = indxt
      go to 10

c     the return is embedded in the above logic

      END

