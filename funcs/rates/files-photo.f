c
c.........................................................................
c

      SUBROUTINE readUVBspec

c     this routine loops through the list of input files and stores the
c     Jnu data

c     Jnu (erg cm^-2 s^-1 Hz^-1 Str^-1) is read in 

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      integer            iz,i
      character*80       infile(UVBzpix),uvbfile,header
      character*255      tabfile      
      include           'rates.com'
      include           'com-photo.com'

      DATA (infile(i),i=1,UVBzpix) / 
     &  'HM05_Jnu_z0.0.cont','HM05_Jnu_z0.2.cont',
     &  'HM05_Jnu_z0.4.cont','HM05_Jnu_z0.6.cont',
     &  'HM05_Jnu_z0.8.cont','HM05_Jnu_z1.0.cont',
     &  'HM05_Jnu_z1.2.cont','HM05_Jnu_z1.4.cont',
     &  'HM05_Jnu_z1.6.cont','HM05_Jnu_z1.8.cont',
     &  'HM05_Jnu_z2.0.cont','HM05_Jnu_z2.2.cont',
     &  'HM05_Jnu_z2.4.cont','HM05_Jnu_z2.6.cont',
     &  'HM05_Jnu_z2.8.cont','HM05_Jnu_z3.0.cont',
     &  'HM05_Jnu_z3.2.cont','HM05_Jnu_z3.4.cont',
     &  'HM05_Jnu_z3.6.cont','HM05_Jnu_z3.8.cont',
     &  'HM05_Jnu_z4.0.cont','HM05_Jnu_z4.2.cont',
     &  'HM05_Jnu_z4.4.cont','HM05_Jnu_z4.6.cont',
     &  'HM05_Jnu_z4.8.cont','HM05_Jnu_z5.0.cont'/
 
      DATA (zuvb(i),i=1,UVBzpix) / 0.0d0,
     &  0.2d0,0.4d0,0.6d0,0.8d0,1.0d0,
     &  1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,
     &  2.2d0,2.4d0,2.6d0,2.8d0,3.0d0,
     &  3.2d0,3.4d0,3.6d0,3.8d0,4.0d0,
     &  4.2d0,4.4d0,4.6d0,4.8d0,5.0d0/
 
c     loop over the UVB data and read in the logJNUuvb matrix

      DO 09 iz=1,UVBzpix

c     construct file name
       uvbfile = infile(iz)
       tabfile = UVBpath
       CALL sappend(tabfile,uvbfile,tabfile)

       OPEN(unit=1,file=tabfile,ERR=999,status='old')
       READ(1,*) header
       READ(1,*) header
       DO 11 i=1,UVBmaxpix
        READ(1,*) logEuvb(i),logJNUuvb(iz,i)
 11    CONTINUE
 
 12    CLOSE(unit=1)

c     set the number of energy data points in the logJNUuvb(iz,i) matrix

 09   CONTINUE

c     set the physical array size 

      NEuvb   = UVBmaxpix

      RETURN

 999  WRITE(6,*) 'ERROR(readUVBspec): UVB spectrum file not found'
      WRITE(6,*) uvbfile
      STOP

      END


c
c.........................................................................
c

      SUBROUTINE readSB99spec

c     this routine loops through the list of input files and stores the
c     Jnu data

c     L_lambda (erg s^-1 Ang^-1) is read in 

c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
 
      include           'rates.h'
      integer            iage,iZ,i
      double precision   lambda,dum
      character*80       infile(SB99Nage,SB99NZ),Sb99file,header  
      character*255      tabfile    
      include           'rates.com'
      include           'com-photo.com'

      DATA (infile(i,1),i=1,SB99Nage) / 
     &  'Ms3_0d004m_1Myr.stb99',
     &  'Ms3_0d004m_5Myr.stb99',
     &  'Ms3_0d004m_10Myr.stb99',
     &  'Ms3_0d004m_20Myr.stb99',
     &  'Ms3_0d004m_40Myr.stb99'/
 
      DATA (infile(i,2),i=1,SB99Nage) / 
     &  'Ms3_1Myr.stb99',
     &  'Ms3_5Myr.stb99',
     &  'Ms3_10Myr.stb99',
     &  'Ms3_20Myr.stb99',
     &  'Ms3_40Myr.stb99'/
 
      DATA (Zsolsb99(i),i=1,SB99NZ) /0.004,1.000/
 
      DATA (Agesb99(i),i=1,SB99Nage) /1.0,5.0,10.0,20.0,40.0/


c     loop over the SB99 data and read in the logLlamsb99 cube


      DO 09 iage=1,SB99Nage

       DO 11 iZ=1,SB99NZ

        Sb99file = infile(iage,iZ)
        tabfile  = Sb99path
        CALL sappend(tabfile,Sb99file,tabfile)

c     read in the file and convert wavelength to energy in eV (we do
c     this only for the first file, because all files have identical
c     wavelengths)

       OPEN(unit=1,file=tabfile,ERR=999,status='old')
       READ(1,'(a80)') header
       READ(1,'(a80)') header
       READ(1,'(a80)') header
       READ(1,'(a80)') header
       READ(1,'(a80)') header
       READ(1,'(a80)') header
       DO 11 i=1,SB99maxpix
        READ(1,*) dum,lambda,logLlamsb99(iage,iZ,i)
        IF (iage.eq.1.AND.iZ.eq.1) then
         logEsb99(i) = log10(h*clight) - log10(lambda) + 8.0d0
     &                 - log10(erg2eV) 
        END IF
 11    CONTINUE
 
 12    CLOSE(unit=1)

 09   CONTINUE

c     set the physical array size

      NEsb99  = SB99maxpix

      RETURN

 999  WRITE(6,*) 'ERROR(readSB99spec): SB99 spectrum file not found'
      WRITE(6,*) SB99file
      STOP

      END


