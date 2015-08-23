
      implicit none

c  integer parameters

c  maxlines   --- max number of lines in the line lists
c  maxpix     --- max number of pixels in a given spectrum
c  maxtrans   --- max number of transitions to make spectra for
c  maxvec     --- max number of pixels for concolution


      integer            maxpix,maxlines,maxtrans,maxvec,maxions
      parameter          (maxlines =  16000      ,  
     @                    maxpix   =  40000      ,
     @                    maxions  =  100        ,
     @                    maxtrans =  1000       ,
     @                    maxvec   =  5.0*maxpix )

      integer            SCREEN,STDIN,STDOUT
      parameter          (SCREEN = 6, STDIN = 5, STDOUT = SCREEN)



c  atomic constants and general acounting

      integer            lines,ndata
      double precision   con1,con2,lambda0,wave_max,wave_min,dwave
      double precision   zabs,snr,vel_min,vel_max,rn,lam0,f0,gamma0
      character*80       ion_name

      COMMON/atomblck/con1(maxlines),con2(maxlines),lambda0(maxlines),
     &            wave_max,wave_min,dwave,zabs,snr,vel_min,vel_max,
     &            rn,lam0(maxlines),f0(maxlines),gamma0(maxlines)
      COMMON/atomchar/ion_name(maxlines)
      COMMON/accntint/lines,ndata



c  data arrays and parameter arrays

      double precision   lambda,wrkflx,sigma
      double precision   zline,nline,bline

      COMMON/fitblck/lambda(maxpix),wrkflx(maxpix),sigma(maxpix)
      COMMON/lineblck/zline(maxlines),nline(maxlines),bline(maxlines)

c  paramlist quantities
     
      integer            nions,ionstage,elemid
      double precision   snrat,vmax
      character*80       element,instr

      COMMON/ipar/nions,ionstage(maxions),elemid(maxions)
      COMMON/rpar/snrat(maxions),vmax(maxions) 
      COMMON/cpar/element(maxions),instr(maxions)


c  convolution arrays 
c  ***** IMPORTANT NOTE FOR CONVOLUTION *****
c  NOTE:  there is a DATA statement in instruments.f that contains
c  a list of powers of 2.  It goes up to MAXCON.  
c  Parameter MAXCON must be a power of 2.

      logical           convolving
      integer           maxcon
      parameter         (MAXCON = 1048576)
      integer           nresponse,ncondat,nfft
      double precision  response,convdata
      double precision  R_fac,slit,conwindo,resfac,profile,hdv,
     &                  pixpres,dv

      COMMON/convlblck/convolving
      COMMON/convrblck/response(maxpix),convdata(maxcon)
      COMMON/conviblck/nresponse,ncondat,nfft
      COMMON/instblck/R_fac,slit,conwindo,resfac,profile,hdv,
     &                pixpres,dv
c
c
c...................................... end minfit.h .................
