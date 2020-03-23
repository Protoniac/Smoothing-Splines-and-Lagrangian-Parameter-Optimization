      SUBROUTINE SPLINES(Y,LDY,H,LDH,Q,LDQ,T,LDT,
     $                   SIG,LDSIG,N,P,A,LDA,B,LDB,C,LDC,
     $                   D,LDD)
!variables données
      DOUBLE PRECISION Y(LDY),H(LDH)
      DOUBLE PRECISION Q(LDQ,*), T(LDT,*), SIG(LDSIG,*)
      INTEGER N
      DOUBLE PRECISION P
!variables résultats 
      DOUBLE PRECISION A(LDA), B(LDB), C(LDC), D(LDD)
!variables locales 
      DOUBLE PRECISION QTQ(N-1,N-1),QSIG(N-1,N+1),M(N-1,N-1)
      DOUBLE PRECISION L(N-1,N-1), TEMP(N-1), SIGQ(N+1,N-1)
!début 
      CALL DGEMM('T','N',N-1,N+1,LDQ,1D0,Q,LDQ,SIG,LDSIG,
     $           0D0,QSIG,N-1)
      CALL DGEMM('N','N',N-1,N-1,LDQ,1D0,QSIG,N-1,Q,LDQ,
     $           0D0,QTQ,N-1)
      CALL ADDITION_MATRICE_CARRE(QTQ,N-1,T,LDT,M,N-1,P)
      CALL CHOLESKY(M,N-1,L,N-1)
      CALL DGEMV('T',LDQ,N-1,P,Q,LDQ,Y,1,
     $           0D0,TEMP,1)
      CALL DESCENTE(L,N-1,TEMP,N-1)
      CALL REMONTEE(L,N-1,TEMP,N-1)
      CALL DGEMM('N','N',LDSIG,N-1,LDQ,1D0,SIG,LDSIG,Q,LDQ,
     $           0D0,SIGQ,N+1)
      CALL DGEMV('N',LDQ,N-1,(1/P)*1D0,SIGQ,N+1,TEMP,1,0D0,A,1)
*f2py intent(inplace) a,b,c,d
      DO I=1,LDA
        A(I) = Y(I) - A(I)
      END DO
      C(1) = 0
      C(LDC) = 0
      DO I=1,LDC-2
        C(I+1) = TEMP(I)
      END DO
      DO I=1,LDD
        D(I) = (C(I+1) - C(I))/(3.0*H(I))
      END DO
      DO I=1,LDB
        B(I) = (A(I+1) - A(I))/H(I) - C(I)*H(I) - D(I)*H(I)*H(I)
      END DO
      DO I=1,LDB
      END DO
      
      END SUBROUTINE
