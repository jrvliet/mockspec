c
c.........................................................................
c

      SUBROUTINE commconfig(GZinfile,z,toler,kdo,jdo,GZoutfile)

c     communicate convergence criterion

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      integer           i,k,j,kk,nek
      integer           kdo(Imax),jdo(Imax)
      double precision  z,toler
      character*80      GZinfile,GZoutfile(Imax)  ! gas+metallicity files
      include           'rates.com'
      include           'com-modes.com'



c     communicate cell properties

      WRITE(6,*) ' ' 
      WRITE(6,*) ' PROCESSING INFO'
      WRITE(6,'(1x,a13,2x,a40)')  ' Gas Box  =',GZinfile
      WRITE(6,'(1x,a13,2x,f11.9)') ' redshift =',z
      WRITE(6,'(1x,a20,i3)')      ' No. of iongrids =',noutfiles
      DO 01 i=1,noutfiles
       WRITE(6,'(1x,a12,2x,a50)')  '  iongrid =',GZoutfile(i)
 01   CONTINUE

      WRITE(4,*) ' ' 
      WRITE(4,*) ' PROCESSING INFO'
      WRITE(4,'(1x,a13,2x,a40)')  ' Gas Box  =',GZinfile
      WRITE(4,'(1x,a13,2x,f11.9)') ' redshift =',z
      WRITE(4,'(1x,a20,i3)')      ' No. of iongrids =',noutfiles
      DO 03 i=1,noutfiles
       WRITE(4,'(1x,a12,2x,a50)')  '  iongrid =',GZoutfile(i)
 03   CONTINUE

c     communicate convergance and checical mixture

      WRITE(6,*) ' ' 
      WRITE(6,*) ' SOLUTION CONVERGENCE'
      WRITE(6,'(1x,a13,1pe11.3)') ' tolerance =',toler
      WRITE(4,*) ' ' 
      WRITE(4,*) ' SOLUTION CONVERGENCE'
      WRITE(4,'(1x,a13,1pe11.3)') ' tolerance =',toler

      WRITE(6,*) ' ' 
      WRITE(6,*) ' INCLUDED ATOMIC SPECIES'
      WRITE(6,'(1x,a13,i4)')      ' Nspecies  =',Nspecies
      WRITE(4,*) ' ' 
      WRITE(4,*) ' INCLUDED ATOMIC SPECIES'
      WRITE(4,'(1x,a13,i4)')      ' Nspecies  =',Nspecies

      WRITE(6,600) 'Z','Elem','Levels'
      WRITE(4,600) 'Z','Elem','Levels'
      DO kk=1,Nspecies
       k = kidx(kk)
       WRITE(6,601) k,specID(k),k+1
       WRITE(4,601) k,specID(k),k+1
      ENDDO
     
      RETURN

 600  FORMAT(1x,3x,a3,a8,2x,a6)
 601  FORMAT(1x,3x,i3,5x,a8,i3)

      END


