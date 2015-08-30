c     FITTING PARAMETERS FOR COLLISIONAL IONIZATION RATE COEFFICIENTS

c     the fitting parameters are from Arnaud + Rothenflung (1985,A&AS,60,425)
c     each parameter corresponds to shell s of ion k,j; we also store the
c     number of shells for ion k,j, and the threshold iomnization energies
c     for shell s of ion k.j

c     shmax(k,j) = number of shells of ion k,j
c     IE(k,j,s)  = ionization threshold energy for shell s of ion k,j
c     A(k,j,s)   = AR85 A coefficent for k,j,s
c     B(k,j,s)   = AR85 B coefficent for k,j,s
c     C(k,j,s)   = AR85 C coefficent for k,j,s
c     D(k,j,s)   = AR85 D coefficent for k,j,s

      integer           shmax
      COMMON/icfit/     shmax(Imax,Imax)

      double precision  IE,A,B,C,D
      COMMON/rcfit/     IE(Imax,Imax,Smax),A(Imax,Imax,Smax),
     &                  B(Imax,Imax,Smax),C(Imax,Imax,Smax),
     &                  D(Imax,Imax,Smax)

