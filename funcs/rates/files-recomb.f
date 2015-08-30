c
c.............................................................................
c
      SUBROUTINE readPLtab

c     this routine reads in the modified version of the pl.dat file
c     obtained from Dima Verner's http://www.pa.uky.edu/~verner/rec.html
c     as published in Arnaud+ (1985,A&AS,60,425)

c     the array architecure is that the 1st indice is the atomic number
c     of the species, A, and the 2nd indice, STG, is the ionization
c     stage, where 1 is neutral (i.e. HeI) and 2 is singly ionized
c     (i.e., HeII); the number of ionization stages tabulated in this
c     file is A-1 for all A

c     A is in column 1 of the file
c     STG is in column 6 of the file

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none     
      include           'rates.h'
      integer           i,j,k,A,Ne,stg,lastA
      double precision  Acoef,Bcoef
      character*80      header,species,specion
      character*255     tabfile
      include           'rates.com'
      include           'com-recomb.com'


c     open the file and read them in; the matrices are sparce in that
c     fitting parameters are stored in the array elements A,stg; all
c     other matrix elements are null

      tabfile = tabpath
      CALL sappend(tabfile,'tab-recomb-PL.dat',tabfile)

      OPEN(unit=1,file=tabfile,status='old')
      READ(1,*) header

      DO 11 i=1,Imax*STGmax
       READ(1,*,END=14) A,Ne,Acoef,Bcoef,species,stg,specion
       A_pl(A,stg) = Acoef
       B_pl(A,stg) = Bcoef
 11   CONTINUE

 14   CLOSE(unit=1)


      RETURN

      END


c
c.............................................................................
c
      SUBROUTINE readLDRtab

c     this routine reads in the modified version of the ldr.dat file
c     obtained from  Dima Verner's http://www.pa.uky.edu/~verner/rec.html
c     as published in Nussbaumer+ (1983,A&A,126,75)

c     the array architecure is that the 1st indice is the atomic number
c     of the species, A; and the 2nd indice, STG, is the ionization
c     stage, where 1 is neutral (i.e. HeI) and 2 is singly ionized
c     (i.e., HeII); the 3rd indice, nTreg, designates whether there are
c     two temperature regimes or one

c     A is in column 1 of the file
c     STG is in column 8 of the file
c     nTreg is in column 10

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none     
      include           'rates.h'
      integer           i,j,k,A,Ne,stg,nTr
      double precision  Acoef,Bcoef,Ccoef,Dcoef,Fcoef
      character*80      header,species,specion
      character*255     tabfile
      include           'rates.com'
      include           'com-recomb.com'


c     open the file and read them in; the matrices are sparce in that
c     fitting parameters are stored in the array elements A,stg; all
c     other matrix elements are null

      tabfile = tabpath
      CALL sappend(tabfile,'tab-recomb-LDR.dat',tabfile)

      OPEN(unit=1,file=tabfile,status='old')
      READ(1,*) header

      DO 11 i=1,Imax*STGmax
       READ(1,*,END=14) A,Ne,Acoef,Bcoef,Ccoef,Dcoef,Fcoef,
     &                  species,stg,nTr,specion
       a_ldr(A,stg,nTr) = Acoef
       b_ldr(A,stg,nTr) = Bcoef
       c_ldr(A,stg,nTr) = Ccoef
       d_ldr(A,stg,nTr) = Dcoef
       f_ldr(A,stg,nTr) = Fcoef
       n1stg(A)         = min(stg,n1stg(A))
       nNstg(A)         = max(stg,nNstg(A))
       nTreg(A,stg)     = max(nTreg(A,stg),nTr)
 11   CONTINUE

 14   CLOSE(unit=1)


c     for good measure, zero the starting ionization stage if there are
c     no ionization stages for this partitular element
      
      DO 15 A=1,Imax
       IF (nNstg(A).eq.0) n1stg(A) = 0
 15   CONTINUE
       
      RETURN

      END


c
c.............................................................................
c
      SUBROUTINE readHDRtab

c     this routine reads in the modified version of the hdr.dat file
c     obtained from  Dima Verner's http://www.pa.uky.edu/~verner/rec.html
c     as published in Arnuad+ (1985, A&AS, 60, 425).

c     the array architecure is that the 1st indice is the atomic number
c     of the species, A, and the 2nd indice, STG, is the ionization
c     stage, where 1 is neutral (i.e. HeI) and 2 is singly ionized
c     (i.e., HeII); the number of ionization stages tabulated in this
c     file is A-1 for all A

c     A is in column 1 of the file
c     STG is in column 6 of the file

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none     
      include           'rates.h'
      integer           i,j,k,A,Ne,stg,lastA
      double precision  Acoef,Bcoef,T0coef,T1coef
      character*80      header,species,specion
      character*255     tabfile
      include           'rates.com'
      include           'com-recomb.com'


c     open the file and read them in; the matrices are sparce in that
c     fitting parameters are stored in the array elements A,stg; all
c     other matrix elements are null

      tabfile = tabpath
      CALL sappend(tabfile,'tab-recomb-HDR.dat',tabfile)

      OPEN(unit=1,file=tabfile,status='old')
      READ(1,*) header

      DO 11 i=1,Imax*STGmax

       READ(1,*,END=14) A,Ne,Acoef,Bcoef,T0coef,T1coef,
     &                  species,stg,specion
       A_hdr(A,stg)  = Acoef
       B_hdr(A,stg)  = Bcoef
       T0_hdr(A,stg) = T0coef
       T1_hdr(A,stg) = T1coef
 11   CONTINUE

 14   CLOSE(unit=1)

      RETURN
      END


