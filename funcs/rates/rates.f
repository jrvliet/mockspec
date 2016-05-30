c.........................................................................
c

      PROGRAM rates

c     this code reads in hydroART data cubes (output from ANA) and
c     computes the ionization conditions in each gas cell as a
c     post-processing step

c     the processing is controlled by the user defined data contained in
c     the file "rates.inp" and "rates.outfiles"; a record of the
c     processing is streamed to screen and written to a file called
c     "rates.runlog"

c     the output is/are the files listed by the use in "rates.outfiles"
c     the output files are data cubes in the same format as the input
c     data cube for each of the ions (referenced k,j, where k is the
c     species index and j is the ionization stage)

c     much of the code physics is described in the paper Churchill etal
c     (arXiv:1409:0916).  

c     Cuurent Version: 3.0
c     Versions 
c     V1.0 optically thin approximation for all cells
c     V2.0 code streamlined for faster computations
c     V3.0 includes self-shielding for optically thick cells

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      include           'getcube.h'
      logical           opticallythick,oflag
      integer           Nthick,Nthin
      integer           ilev,i,j,k,kk
      integer           update,report
      integer           kdo(Imax),jdo(Imax)
      double precision  toler
      double precision  z,T,Zmet
      character*80      GZinfile,GZcoolfile,GZoutfile(Imax) 
      include           'rates.com'
      include           'com-modes.com'
      include           'getcube.com'
      character*300     codeLoc
      character*20      tab,uvb,sb




      OPEN(unit=4,file='rates.runlog',status='unknown')

      CALL notify(' ',-1)
      CALL notify('........................................',-1)
      CALL notify('.                                      .',-1)
      CALL notify('.         GRID ION CODE V3.0           .',-1)
      CALL notify('.         with selfshielding           .',-1)
      CALL notify('.           November 2014              .',-1)
      CALL notify('.                                      .',-1)
      CALL notify('........................................',-1)
      CALL notify(' ',-1)


      i       = 0
      Nprint  = 0
      oflag   = .false.

c     table, UVB, and Sb99 data paths (global)
c     TABPATH is where the atomic data tables live
c     UVBPATH is where the Haardt & Madau UVB SEDs live
c     SB99PATH is where the Starburst99 SEDs live
c     Get this information from the command line
      CALL getarg(1, codeLoc)

c     Build the path to the tables
      tab = '/data/grid/'
      uvb = '/data/uvb/'
      sb = '/data/sb99/'
c      tabpath  = trim(codeLoc) // trim(tab)
c      UVBpath  = trim(codeLoc) // trim(uvb)
c      Sb99path = trim(codeLoc) // trim(sb)

c      tabpath='/lustre/projects/p089_swin/jvander/mockspec/data/grid/'
c      UVBpath='/lustre/projects/p089_swin/jvander/mockspec/data/uvb/'
c      Sb99path='/lustre/projects/p089_swin/jvander/mockspec/data/sb99'

c      tabpath  = '/home/matrix2/cwc/Projects/Mockspec/GridCode/'
c      UVBpath  = '/home/matrix2/cwc/Projects/Mockspec/UVBspectrum/'
c      Sb99path = '/home/matrix2/cwc/Projects/Mockspec/Sb99spectrum/'

c      tabpath  = '/home/jacob/research/code/mockspec/data/grid/'
c      UVBpath  = '/home/jacob/research/code/mockspec/data/uvb/'
c      Sb99path = '/home/jacob/research/code/mockspec/data/sb99/'

      tabpath  = '/home/hyades/jrvander/mockspec/data/grid/'
      UVBpath  = '/home/hyades/jrvander/mockspec/data/uvb/'
      Sb99path = '/home/hyades/jrvander/mockspec/data/sb99/'

c     ............................................
c     PHASE 1: SET UP THE RUN, COMMUNICATE TO USER
c     ............................................

c     1. configure arrays, globals, zero out all arrays, etc
c     2. read in the "rates.outfiles" file and configure the arrays that
c        track which ions are to be output for processing
c     3. communicate the current run to the user
c     4. configure the physics
c     5. read in the data cube
c     6. open the output cubes for the target ions

      CALL configall(toler,GZinfile,z)
      CALL outfiles(kdo,jdo,GZoutfile,GZcoolfile)
      CALL commconfig(GZinfile,z,toler,kdo,jdo,GZoutfile)
      CALL configIonPhysics(z)
      CALL readGZ(GZinfile)
      CALL openfiles(GZoutfile,GZcoolfile)

c     ...................................
c     PHASE 2: POST-PROCESS THE DATA CUBE
c     ...................................

c     make the SED (this will need to be moved in the logical flow of
c     the program when we incorporate stars); and compute the
c     photoionization rates for each shell, then populate the
c     photoionization and Auger rates

      CALL mkSED
      CALL save_photoRs

c     loop over cell levels: if doing self shielding create the
c     ionization grid; if all cells in this level are optically thin
c     then OFLAG remains set low and we solve the cells directly; if
c     OFLAG is returned high, we employ self-shielding (note that if we
c     are not doing self-shielding OFLAG will remain in a low state and
c     solvethick will never be called)

      DO ilev=1,Nclevels

       CALL notify(' SOLVING CELLS: Level =',ilev)

       Nthin  = 0
       Nthick = 0

       IF (doPH.AND.doSLFSH) then ! engage self-shielding
        CALL Rkjsgrid(ilev,kdo,jdo,toler,oflag)
        IF (oflag) then ! apply self-shielding as needed
         CALL solvethick(ilev,kdo,jdo,Nthin,Nthick,z,toler)  
        ELSE ! self-shielding not needed for this level
         CALL solvethin(ilev,kdo,jdo,Nthin,z,toler) 
        ENDIF
       ELSE ! not doing shielding
        CALL solvethin(ilev,kdo,jdo,Nthin,z,toler) 
       ENDIF

       CALL notify('   No. of optically thin cells =',Nthin)
       CALL notify('   No. of self-shielded cells  =',Nthick)

      ENDDO ! next level

c     close all output files

      CALL notify('  ALL CELLS PROCESSED',-1)
      CALL notify('  closing files',-1)
      CALL closefiles
      
c     communicate completion

      CALL notify(' ',-1)
      CALL notify(' SUCCESSFUL COMPLETION',-1)
    
c     ....
c     DONE
c     ....

c     close the runlog file

      CLOSE(unit=4)

      STOP
      END

c
c.........................................................................
c

      SUBROUTINE solvethin(ilev,kdo,jdo,Nthin,z,toler) 

c     if all the cells on cell level ILEV are optically thin, then we
c     solve them directly

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
  
      implicit none
      include           'rates.h'
      include           'getcube.h'
      integer           Nthin
      integer           ilev,i,m
      integer           k,j,kk
      integer           kdo(Imax),jdo(Imax)
      double precision  z,T,Zmet
      double precision  toler
      include           'rates.com'
      include           'com-modes.com'
      include           'getcube.com'
          

c     communicate

      CALL notify('  number of cells = ',Nclevel(ilev))
      CALL notify('   ID start       = ',ID1(ilev))
      CALL notify('   ID end         = ',IDN(ilev))

c     for good measure, populate the photoionization and Auger rates

      CALL pop_photoRs

c     loop over the cells for this cell level

      DO 01 i=ID1(ilev),IDN(ilev)

c     for cell i: store nH, T, and Zmet

       nH   = den(i)
       T    = Tcell(i)
       Zmet = SNII(i) + SNIa(i)

c     compute abundance pattern
c     compute the collisional rate coefficients
c     solve ionization balance
c     compute the ionization and recombination timescales

      Nthin = Nthin + 1
      CALL initabund(Zmet)
      CALL alphas_T(T)     
      CALL ionbalance(toler)
      CALL timescales(T,Zmet,kdo,jdo)

c     write the cell data to the new ion cubes

      DO 92 m=1,noutfiles
       CALL writeGZcell(i,m,kdo,jdo)
 92   CONTINUE

c     write the data required to compute the photo cooling times, which
c     is done outside this program post processing

      CALL writecooldat(i,z,T,Zmet) 

c     go to the next cell i

 01   CONTINUE

      RETURN
      END

c
c.........................................................................
c

      SUBROUTINE solvethick(ilev,kdo,jdo,Nthin,Nthick,z,toler)  

c     if some of the cells on cell level ILEV are optically thick, then
c     we solve the optically thick cells using the grid, but when we
c     encounter an optically thin cells, we solve them directly

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      include           'rates.h'
      include           'getcube.h'
      logical           tempflag
      integer           Nthin,Nthick
      integer           ilev,i,m
      integer           k,j,kk
      integer           kdo(Imax),jdo(Imax)
      double precision  z,T,Zmet
      double precision  toler
      include           'rates.com'
      include           'com-modes.com'
      include           'getcube.com'
          

c     communicate

      CALL notify('  number of cells = ',Nclevel(ilev))
      CALL notify('   ID start       = ',ID1(ilev))
      CALL notify('   ID end         = ',IDN(ilev))

c     for good measure, populate the photoionization and Auger rates

      CALL pop_photoRs

c     loop over the cells for this cell level

      DO 01 i=ID1(ilev),IDN(ilev)

c     for cell i: store nH, T, and Zmet

       nH   = den(i)
       T    = Tcell(i)
       Zmet = SNII(i) + SNIa(i)

c     compute abundance pattern

       CALL initabund(Zmet)

c     determine if the cell requires self shielding treatment or is
c     optically thin; we determine this by checking if the gas phase
c     location for this cell level resides in the regime of optically
c     thick cells for this level

       IF    ((log10(T).ge.Tmin(ilev))
     &  .AND. (log10(T).le.Tmax(ilev))
     &  .AND. (log10(nH).ge.nHmin(ilev))
     &  .AND. (log10(nH).le.nHmax(ilev)) ) then

        tempflag = .true.   ! TEMPORARY
        Nthick = Nthick + 1
        CALL selfshield(T,kdo,jdo)

C**     In order to do the timescales, we either need to make a grid of
C**     the three timescales, or a grid of the photoionization rate, and
C**     the recomination and collisional ionization rate coeffs - boo

       ELSE

c     if the cell is optically thin, then solve it for the unattenuated
c     SED; compute the collisional rate coefficients solve ionization
c     balance

        tempflag = .false.  ! TEMPORARY
        Nthin = Nthin + 1
        CALL alphas_T(T)     
        CALL ionbalance(toler)
   
       END IF

c     TEMPORARY: treatment of timescales here
c     write the cell data to the new ion cubes

       CALL timescales(T,Zmet,kdo,jdo)  ! cooling time only (NEEDS FIX)
        DO 92 m=1,noutfiles
         IF (tempflag) then   ! TEMPORARY 
          tau_ph(m)   = 0.0d0
          tau_rec(m)  = 0.0d0
          tau_coll(m) = 0.0d0
         ENDIF
        CALL writeGZcell(i,m,kdo,jdo)
 92    CONTINUE

c     write the data required to compute the photo cooling times, which
c     is done outside this program post processing

       CALL writecooldat(i,z,T,Zmet) 

c     go to the next cell i

 01   CONTINUE

      RETURN
      END


c
c.........................................................................
c

c     INCLUDE MODULES (subroutines and functions)

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

c     data cube I/O routines

      include 'getcube.f'

c     communication routines and general file I/O

      include 'rates-comm.f'
      include 'rates-files.f'

c     computation of SED, rate coefficients, photo cross sections

      include 'rates-sed.f'
      include 'rates-alphas.f'
      include 'rates-phxsecs.f'
      include 'rates-Rkjs.f'

c     photoionization calculations and population of tabulated data

      include 'funcs-photo.f'   
      include 'files-photo.f'
      include 'phfit-mod.f'
      include 'phfit-data.blk'

c     recombination calculations and population of tabulated data

      include 'funcs-recomb.f'
      include 'files-recomb.f'

c     collisional calculations and population of tabulated data

      include 'funcs-coll.f'
      include 'files-coll.f'

c     charge transfer calculations and population of tabulated data

      include 'kingdon-mod.f'
      include 'kingdon-data.blk'

c     auger calculations and population of tabulated data

      include 'files-auger.f'
      include 'funcs-auger.f'

c     general subs, rate matrix solver, gas timescales

      include 'rates-subs.f'
      include 'rates-solve4.f'
      include 'rates-thick.f'
      include 'rates-times.f'

c     eof
