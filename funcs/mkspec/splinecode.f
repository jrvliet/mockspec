c     the Numerical Recipes splining routines used for the
c     interpolation of the ionization correction tables
c
c     this file contains:
c     SUBROUTINE SPLIE2
c     SUBROUTINE SPLIN2
c     SUBROUTINE SPLINE
c     SUBROUTINE SPLINT
c
c
c     ROUTINE SPLIE2 computes 2D 2nd derivatives of grid

      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=10000)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN)

      DO 13 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
        CALL SPLINE(X2A,YTMP,N,1.D30,1.D30,Y2TMP)
        DO 12 K=1,N
          Y2A(J,K)=Y2TMP(K)
12      CONTINUE
13    CONTINUE

      RETURN
      END

c     ROUTINE SLIN2 performs 2D interpolation

      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=10000)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),
     &          Y2TMP(NN),YYTMP(NN)

      DO 12 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
          Y2TMP(K)=Y2A(J,K)
11      CONTINUE
        CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J))
12    CONTINUE

      CALL SPLINE(X1A,YYTMP,M,1.D30,1.D30,Y2TMP)
      CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y)

      RETURN
      END

c     ROUTINE SPLINE computes 1D 2nd derivatives of grid

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=10000)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)

      IF (YP1.GT..99D30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF

      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE

      IF (YPN.GT..99D30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)

      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE

      RETURN
      END

c     ROUTINE SPLINT performs 1D interpolation

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(N),YA(N),Y2A(N)

      KLO=1
      KHI=N

1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

      H=XA(KHI)-XA(KLO)

      IF (H.EQ.0.) then
       STOP 'Bad XA input in Num Recipe routine SPLINT.'
      END IF

      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.

      RETURN
      END
