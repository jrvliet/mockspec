
      integer            nmx,mxlin,mxord,mxions,mxtrans
      parameter          (nmx     =  10000, 
     &                    mxlin   =     20, 
     &                    mxord   =     15,
     &                    mxions  =     30,
     &                    mxtrans =   300)

      integer            M,Jnot
      parameter          (Jnot = 6, M = 2*Jnot + 1)

      double precision   toler
      parameter         (toler = 1.0d-6)
