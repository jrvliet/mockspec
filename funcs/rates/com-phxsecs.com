c     PHOTOIONIZATION PARTIAL CROSS SECTION DATA

c     global variables for photoionization compuations; the cross sections
c     are a function of photon energy

c     NEsigkjs    = number of array elements in each k,j,s array 
c     logEsigkjs  = photon energy [eV] array for logsigkjs andd d2sigkjsdE2 
c     logsigkjs   = log10 partial cross sections for shell s of ion k,j
c     d2sigkjsdE2 = 2nd derivitives for spline interpolation of logsigkjs

      integer           NEsigkjs
      COMMON/Exkjs/     NEsigkjs(Imax,Imax,Nshmax)

      double precision  logEsigkjs,logsigkjs,d2sigkjsdE2
      COMMON/xseckjs/   logEsigkjs(NEmax,Imax,Imax,Nshmax),
     &                  logsigkjs(NEmax,Imax,Imax,Nshmax),
     &                  d2sigkjsdE2(NEmax,Imax,Imax,Nshmax)
