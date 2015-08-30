c     FITTING PARAMETERS FOR RECOMBINATION PHYSICS

c     includes two-body radiative recombination (PL), low temperature
c     dielectronic recombination (LDR), and high temperature dielectronic 
c     recombination (HDR)


c     PL DATA
      integer           Ne_pl
      COMMON/iPLtab/    Ne_pl(Imax,STGmax)

      double precision  A_pl,B_pl
      COMMON/rPLtab/    A_pl(Imax,STGmax),
     &                  B_pl(Imax,STGmax)


c     LDR DATA
      integer           n1stg,nNstg,nTreg
      COMMON/iLDRtab/   n1stg(Imax),nNstg(Imax),
     &                  nTreg(Imax,STGmax)

      double precision  a_ldr,b_ldr,c_ldr,d_ldr,f_ldr
      COMMON/rLDRtab/   a_ldr(Imax,STGmax,2),
     &                  b_ldr(Imax,STGmax,2),
     &                  c_ldr(Imax,STGmax,2),
     &                  d_ldr(Imax,STGmax,2),
     &                  f_ldr(Imax,STGmax,2)


c     HDR DATA
      double precision  A_hdr,B_hdr,T0_hdr,T1_hdr
      COMMON/rHDRtab/   A_hdr(Imax,STGmax),
     &                  B_hdr(Imax,STGmax),
     &                  T0_hdr(Imax,STGmax),
     &                  T1_hdr(Imax,STGmax)

