
      logical            sf_flag,int_flag
      integer            ndata,nfind,f_beg,f_end,mask,norders,
     &                   collim,nlines
      integer            ntie,tie,jtie,njtie
      integer            nions,ionstage
      double precision   fval,lambda0,zbar
      double precision   wave,vel,flux,sigma,smsig,cont,
     &                   tau,dutau,ddtau,col,ducol,ddcol,
     &                   ewpix,ewsigpix
      double precision   ew,ewsig,wbar,sigwbar,width,sigwidth,
     &                   dr,drsig,siglevel
      double precision   vbar,vwidth,vasym,vtotbar,vtotwidth,
     &                   vtotasym,sigvbar,sigvwidth,sigvasym,
     &                   sigvtotbar,sigvtotwidth,sigvtotasym
      double precision   ewtot,ewsigtot,drtot,drsigtot,
     &                   sltot,tautot,dutautot,ddtautot,
     &                   coltot,ducoltot,ddcoltot,N_sigma,
     &                   Nsigion,Slevid
      double precision   Slev,profile
      double precision   EWlim,EWcut
      double precision   wave0,fosc
      character*80       order,ionfile,spec_name
      character*80       specfiles,element,instr,instrid


      COMMON /lllblk/  sf_flag(mxord,mxlin),int_flag

      COMMON /tieblk/  ntie,tie(mxord),jtie(mxord,mxord),
     &                 njtie(mxord)

      COMMON /iiiblk/  ndata(mxord),nfind(mxord,mxlin),
     &                 f_beg(mxlin,mxord,mxlin),
     &                 f_end(mxlin,mxord,mxlin),
     &                 mask(nmx,mxord),norders,
     &                 collim(nmx,mxord),nlines,
     &                 ionstage(mxions),nions

      COMMON /rrdblk/ wave(nmx,mxord),vel(nmx,mxord),
     &                flux(nmx,mxord),sigma(nmx,mxord),
     &                smsig(nmx,mxord),cont(nmx,mxord),
     &                ewpix(nmx,mxord),ewsigpix(nmx,mxord),
     &                tau(nmx,mxord),dutau(nmx,mxord),
     &                ddtau(nmx,mxord),col(nmx,mxord),
     &                ducol(nmx,mxord),ddcol(nmx,mxord),
     &                EWlim(mxions),EWcut(mxions)

      COMMON /rrlblk/ ew(mxlin,mxord),ewsig(mxlin,mxord),
     &                wbar(mxlin,mxord),sigwbar(mxlin,mxord),
     &                width(mxlin,mxord),sigwidth(mxlin,mxord),
     &                dr(mxlin,mxord),drsig(mxlin,mxord),
     &                siglevel(mxlin,mxord),vbar(mxlin,mxord),
     &                vwidth(mxlin,mxord),vasym(mxlin,mxord),
     &                sigvbar(mxlin,mxord),sigvwidth(mxlin,mxord),
     &                sigvasym(mxlin,mxord)


      COMMON /rroblk/ ewtot(mxord),ewsigtot(mxord),
     &                drtot(mxord),drsigtot(mxord),
     &                coltot(mxord),ducoltot(mxord),
     &                ddcoltot(mxord),tautot(mxord),
     &                dutautot(mxord),ddtautot(mxord),
     &                sltot(mxord),N_sigma(mxord),
     &                vtotbar(mxord),vtotwidth(mxord),
     &                vtotasym(mxord),sigvtotbar(mxord),
     &                sigvtotwidth(mxord),sigvtotasym(mxord),
     &                Slev,profile,Nsigion(mxions),Slevid(mxions)

      COMMON /miscblk/ fval(mxord),lambda0(mxord),zbar,
     &                 wave0(mxions,mxtrans),
     &                 fosc(mxions,mxtrans)

      COMMON /cccblk/  order(mxord),ionfile,spec_name(mxord),
     &                 specfiles(mxord),element(mxions),
     &                 instr(mxions),instrid(mxions)
