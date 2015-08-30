c     BOOLEAN PHYSICS MODES

c     these booleans are entered into the rates.inp file by the use
c     when they are set high, then the corresponding ionization physics 
c     is included, when set low, the physics is not included

c     doPH     = photoionization
c     doLDR    = lo T dielectronic recombination
c     doHDR    = hi T dielectronic recombination
c     doCDI    = direct collisional ionization
c     doCEA    = collisional excitation-autoionization 
c     doAUG    = auger ionization (which is coupled to photoionization)
c     doCT     = charge transfer ionizations and recombinations
c     ubvflag  = UVB ionizing spectrum
c     sb99flag = SB99 stellar population ionizing spectrum
c     doSLFSH  = self shielding of optically thick cells

      logical           doPH,doLDR,doHDR,doCDI,doCEA,doAUG,doCT,
     &                  uvbflag,sb99flag,doSLFSH
      COMMON/physmodes/ doPH,doLDR,doHDR,doCDI,doCEA,doAUG,doCT,
     &                  uvbflag,sb99flag,doSLFSH

