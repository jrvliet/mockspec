
      integer           nlos,nsys,nsubsys,nreg,nsubreg
      COMMON/iblock/    nlos,nsys,nsubsys,
     &                  nreg(maxlos),nsubreg(maxlos)

      double precision  impact
      COMMON/rblock/    impact(maxlos)

      character*80      losfile,sysabsfile,regabsfile,losID
      COMMON/cblock/    losfile(maxlos),sysabsfile(maxlos),
     &                  regabsfile(maxlos),losID(maxlos)

      character*250     headline,sysdata,regdata
      COMMON/cblk2/     headline,sysdata(maxlos,mxreg),
     &                  regdata(maxlos,mxreg)
