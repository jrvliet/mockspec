c..............................................................................
c
      PROGRAM             bindata

c     performs log binning of an arbitrary data array
c     Usage: bindata2 $1 $2 $3 $4 $5 $6 $7
c       $1  = the input file name 
c       $2  = column the data are in 
c       $3  = input data type (0=linear,1=log10)
c       $4  = binsize (equal log10)
c       $5  = lower bin limit (log10)
c       $6  = upper bin limit (log10)
c       $7  = number of header rows to skip

c     the data file should be formated as shown, list directed in
c     (row,column) order as follows

c       hdr11   hdr12       hdr1m
c       hdr21   hdr22       hdr2m
c       x(1,1) x(1,2) . . . x(1,m)
c       x(2,1) x(2,2) . . . x(2,m)
c       .  .                .     
c       .  .                .
c       .  .                .
c       x(n,1) x(n,2) . . . x(n,m)
 
c..............................................................................
  
      implicit none
      external          CLup,CLdn
      integer           maxdata
      parameter         (maxdata=400000)
      integer           i,j,col,ndata,nbins,type,nhdr
      double precision  zbrent
      double precision  n,Cl,tol,xllim,xulim,lgxcen,lgxlo,lgxhi
      double precision  x(maxdata),bin(maxdata)
      double precision  binup(maxdata),bindn(maxdata)
      double precision  binsize,z1,z2,lo,hi,xx
      character*400     in_file,out_file,string
      character*400     headerline
 
      COMMON /poisson/  n,CL

c     we assume a 1-sigma single sided confidence level

      CL  = 0.8413
 
c     we assume a tolerance of 1.e-5 for the numerical root solving for
c     the Poisson error; the maximum value of any bin is XULIM; if we
c     exceed this we will need to catch it and terminate

      tol   = 1.0d-7

c     get the command line arguments

      CALL commline(in_file,col,type,binsize,z1,z2,nhdr)


c      WRITE(6,*) 'Input file             ',in_file
c      WRITE(6,*) 'Data column            ',col
c      WRITE(6,*) 'bin width              ',binsize
c      WRITE(6,*) 'start domain value     ',z1
c      WRITE(6,*) 'end domain star value  ',z2
c      WRITE(6,*) 'header rows            ',nhdr

      CALL fappend(in_file,'logfreqbin',out_file)

c      WRITE(6,*) 'OUTPUT FILE NAME = ',out_file


c     enter main program here

  
c     -----------------------------------------------------------------
c     READ in the data file; store data in x(i); count the data (NDATA)
c     -----------------------------------------------------------------

      ndata = 0
c      write(*,*) in_file
      OPEN(unit=1,file=in_file,ERR=998,status='old')
      IF (nhdr.gt.0) then
c       WRITE(6,*) 'skipping ',nhdr,' header rows ...'
       DO 09 i=1,nhdr
        READ(1,*,ERR=999) headerline
 09    CONTINUE
      END IF
      DO 11 i=1,maxdata
       READ(1,*,ERR=999,END=12) (xx,j=1,col-1),x(i)
       If (type.eq.1) x(i) = 10.0d0**x(i)  ! convert to linear for binning   
       ndata = ndata + 1                   ! integer number of data points
 11   CONTINUE
 12   CLOSE(unit=1)
  
c     -----------------------------------------------------------------
c     bin the data; datum is in a bin if it is >= the lower bound and <
c     the upper bound, but don't miss xmax; store normalization
c     possibilities on the fly (MAXBIN,NORM); compute the uncertainty in
c     the bin value (use Poisson stats)
c     -----------------------------------------------------------------

c     set counters

      nbins = nint((z2-z1)/binsize)

c     loop over bins

      DO 15 i=1,nbins 

c     initialize

       n      = 0.0d0
       bin(i) = 0.0d0
       lo     = 10.0d0**z1 * (10.0d0**((i-1)*binsize))
       hi     = 10.0d0**z1 * (10.0d0**(i*binsize))

c     count how many are in bin; null uncertainties 

       DO 13 j=1,ndata
        IF ((x(j).GE.lo).AND.(x(j).LT.hi)) bin(i)=bin(i)+1.0d0
        IF ((i.EQ.nbins).AND.(x(j).EQ.hi)) bin(i)=bin(i)+1.0d0
        bindn(i) = 0.0d0
        binup(i) = 0.0d0
 13    CONTINUE

c     compute Poisson errors

c     if bin is large, then use Gaussian approx

       IF (bin(i).gt.200.0) then                 

        binup(i) = bin(i) + sqrt(bin(i))        
        bindn(i) = bin(i) - sqrt(bin(i))        

c     otherwise use Possion CDF to compute uncertainties, which only can
c     happen if bin(i).ne.0

       ELSE                                    

        If (bin(i).gt.0.0) then
         n         = bin(i)                ! the value to be passed
         xllim     = 0.0d0                 ! minimum domain for root solve
         xulim     = 1.0d3*n               ! maximum domain for root solve
c         WRITE(6,*) i,xllim,xulim
         bindn(i) = zbrent(CLdn,xllim,xulim,tol) ! dn Poisson mean
         binup(i) = zbrent(CLup,xllim,xulim,tol) ! up Poisson mean
        END IF         

       END IF

 15   CONTINUE
  

c     -----------------------------------------------------------------
c     normalize and output
c     -----------------------------------------------------------------

      OPEN(unit=2,file=out_file,status='unknown')
      WRITE(2,199) in_file
      WRITE(2,200) 

      DO 23 i=1,nbins

c     binning linear and log10

       lo     = 10.0d0**z1 * (10.0d0**((i-1)*binsize))
       hi     = 10.0d0**z1 * (10.0d0**(i*binsize))
       lgxlo  = log10(lo)
       lgxhi  = log10(hi)
       lgxcen = 0.50d0*(lgxhi+lgxlo)

c     frequency distribution normalization (per x variable per Ndata)

       bin(i)   = bin(i)/(hi-lo)/float(ndata)   
       bindn(i) = bindn(i)/(hi-lo)/float(ndata)
       binup(i) = binup(i)/(hi-lo)/float(ndata)

c     convert errors to log10

       IF (bin(i).gt.0.0) then
         bindn(i) = log10(bindn(i)/bin(i))    
         binup(i) = log10(binup(i)/bin(i))   
       ELSE
         bindn(i) = 0.0d0
         binup(i) = 0.0d0
       END IF

c     convert X data to log10

       IF (bin(i).gt.0.0) then
        bin(i) = log10(bin(i))
       ELSE
        bin(i) = -50.0
       END IF

       WRITE(2,201) lgxcen,bin(i),0.5*(lgxhi-lgxlo),bindn(i),binup(i)
 23   CONTINUE
      CLOSE(unit=2)

c     fini

      STOP 

c     formats
  
 199  FORMAT(1x,/,1x,'Input data file: ',a)
 200  FORMAT(1x,/,1x,t5,'bin center',t18,'log(F)',
     @       t31,'|bin/2|',t44,'dlogFdn',t57,'dlogFup')
 201  FORMAT(1x,1pe13.5,1p5e13.5)

c     error trapping

 998  WRITE(6,*) '  ************************'
      WRITE(6,*) '  * Input file not found *'
      WRITE(6,*) '  ************************'
      STOP 'ERROR(bindata): input file cannot be opened'
  
 999  WRITE(6,*) '  *************************************' 
      WRITE(6,*) '  * Input file format is inconsistent *'
      WRITE(6,*) '  * READ or re-READ the instructions. *'
      WRITE(6,*) '  *************************************'
      STOP 'ERROR(bindata): bad col number or input file format'
  
      END
  
c.............................................................................
  
      SUBROUTINE          helpline

c.............................................................................

      WRITE(6,*) ' BINDATA-LOGFREQ bins a column of data.  The X data'
      WRITE(6,*) ' can be linear or log- the user specifies.  The '
      WRITE(6,*) ' binning is performed with a constant log10 bin size.'
      WRITE(6,*) ' The binned values are output in as log10 frequency'
      WRITE(6,*) ' which is to say that the counts in each bin are '
      WRITE(6,*) ' normalized by the linear bin width. Uncertainties'
      WRITE(6,*) ' are computed from a Poisson density function using'
      WRITE(6,*) ' a single-sided confidence level of 1 sigma.  Thus,'
      WRITE(6,*) ' error bars on the frequencies may be asymmetric.'
      WRITE(6,*) ' Usage: bindata2 $1 $2 $3 $4 $5 $6 $7'
      WRITE(6,*) '  $1  = the input file name '
      WRITE(6,*) '  $2  = column the data are in' 
      WRITE(6,*) '  $3  = input data type (0=linear,1=log10)'
      WRITE(6,*) '  $4  = binsize (equal log10)'
      WRITE(6,*) '  $5  = lower bin limit (log10)'
      WRITE(6,*) '  $6  = upper bin limit (log10)'
      WRITE(6,*) '  $7  = number of header rows to skip'

      RETURN
      END

c..............................................................................

      SUBROUTINE          fappend(in_file,delim,out_file)   

c  replace the delimeter on infile with delim and output the
c  result in out_file
c  If the infile has no delimeter, then append delim prefixed 
c  with a "." to infile and output the result in out_file

c..............................................................................
  
      implicit none
      integer             i,k,lend      
      character*(*)       in_file,delim,out_file
  
       lend = len(in_file)
       k    = 0
       DO 09 i=1,lend
        k = i
        IF ((in_file(i:i).eq.'.').or.(in_file(i:i).eq.' ')) goto 10
 09   CONTINUE
  
 10   IF (in_file(k:k).eq.'.') out_file= in_file(1:k)//delim
      IF (in_file(k:k).eq.' ') out_file= in_file(1:k-1)//"."//delim

      RETURN
      END
  
c     the functions for the Poisson uncertaities

      DOUBLE PRECISION FUNCTION CLup(x)  ! upward uncertainty
      double precision x,gammp,n,CL
      COMMON /poisson/ n,CL
      CLup = CL - gammp(n+1.0d0,x)
      RETURN
      END

      DOUBLE PRECISION FUNCTION CLdn(x)  ! downward uncertainty
      double precision x,gammq,n,CL
      COMMON /poisson/ n,CL
      CLdn = CL - gammq(n,x)
      RETURN
      END



      SUBROUTINE commline(in_file,col,type,binsize,z1,z2,nhdr)

      implicit none
      logical           error
      integer           col,type,nhdr
      double precision  binsize,z1,z2,value2
      character*400     in_file,string
      
c     input file name

      CALL getarg(1,string)  
      IF (string.eq.' ') then
       CALL helpline
       STOP '(bindata-logfreq): comm-line help'
      ELSE
       in_file = string
      END IF

c    column data are in

      CALL getarg(2,string)
      col = value2(string,error)
      IF (error) then
       WRITE(6,*) ' ERROR: bad $2 on comm-line'
       CALL helpline
       STOP 
      END IF

c    are data in log10 already or linear (0=linear, 1=log10)

      CALL getarg(3,string)
      type = int(value2(string,error))
      IF (error) then
       WRITE(6,*) ' ERROR: bad $3 on comm-line'
       CALL helpline
       STOP 
      END IF

c     binsize log10 (equal space in dex)

      CALL getarg(4,string)
      binsize = value2(string,error)
      IF (error) then
       WRITE(6,*) ' ERROR: bad $4 on comm-line'
       CALL helpline
       STOP 
      END IF

c     lower domain limit log10

      CALL getarg(5,string)
      z1 = value2(string,error)
      IF (error) then
       WRITE(6,*) ' ERROR: bad $5 on comm-line'
       CALL helpline
       STOP 
      END IF

c     upper domain limit log10

      CALL getarg(6,string)
      z2 = value2(string,error)
      IF (error) then
       WRITE(6,*) ' ERROR: bad $6 on comm-line'
       CALL helpline
       STOP 
      END IF

c     header lines

      CALL getarg(7,string)
      IF (string.eq.' ') then
       nhdr = 0
      else
       nhdr = int(value2(string,error))
      END IF

      RETURN
      END

c  eof   
