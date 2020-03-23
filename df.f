      SUBROUTINE DFP(Q,LDQ,T,LDT,SIG,LDSIG,N,RES,LDRES,
     $                             IND,LDIND,A,B,H)
! Cette fonction permet de calculer un vecteur de df en fonction de p
      INTEGER LDQ,LDT,LDSIG,N,I,INDICE,P_LOOP
      DOUBLE PRECISION Q(LDQ,*), T(LDT,*), SIG(LDSIG,*)
      DOUBLE PRECISION RES(LDRES),IND(LDIND)
!variables locales
      DOUBLE PRECISION IDENT(N+1,N+1),M(N-1,N-1),
     $                 QSIG(N-1,N+1),SIGQ(N+1,N-1),
     $                 QQ(N-1,N-1),
     $                 INVERSE(N-1,N-1), TEMP(N+1,N-1)
      DOUBLE PRECISION A,B,H,P,DF
      INDICE = 1
      DF = 0
*f2py intent(inplace) res
*f2py intent(inplace) ind
      CALL IDENTITE(IDENT,N+1)
      CALL DGEMM('T','N',LDQ-2,N+1,LDQ,1D0,Q,LDQ,SIG,LDSIG,
     $              0D0,QSIG,N-1)
      CALL DGEMM('N','N',N-1,N-1,N+1,1D0,QSIG,N-1,Q,LDQ,
     $              0D0,QQ,N-1)
      DO P_LOOP = 1,INT(B*H)
        P = A + (P_LOOP - 1)/(H * 1D0)
        CALL ADDITION_MATRICE_CARRE(QQ,N-1,T,LDT,M,N-1,P)
        CALL INVERSION(M,N-1,INVERSE,N-1,N-1)
        CALL DGEMM('N','N',LDSIG,N-1,N+1,1D0,SIG,LDSIG,Q,LDQ,
     $              0D0,SIGQ,N+1)
        CALL DGEMM('N','N',N+1,N-1,N-1,1D0,SIGQ,N+1,INVERSE,
     $              N-1,0D0,TEMP,N+1)
        CALL DGEMM('N','T',N+1,N+1,N-1,1D0,TEMP,N+1,Q,
     $              LDQ,0D0,IDENT,N+1)
        DO I = 1,N+1
            DF = DF + (1-IDENT(I,I))
        END DO
        RES(INDICE) = DF
        IND(INDICE) = P
        INDICE = INDICE + 1
        DF = 0D0
      END DO
      END SUBROUTINE
      
      DOUBLE PRECISION FUNCTION PDF(Q,LDQ,T,LDT,SIG,
     $                              LDSIG,N,DFGOAL,A,B)
! Cette fonction permet de calculer p selon un df donné par
! recherche dichotomique
      INTEGER LDQ,LDT,LDSIG,N,I,INDICE
      DOUBLE PRECISION Q(LDQ,*), T(LDT,*), SIG(LDSIG,*)
      DOUBLE PRECISION SEUIL,DF,A,B,DFGOAL
!variables locales
      DOUBLE PRECISION IDENT(N+1,N+1),M(N-1,N-1),
     $                 QSIG(N-1,N+1),SIGQ(N+1,N-1),
     $                 QQ(N-1,N-1),
     $                 INVERSE(N-1,N-1), TEMP(N+1,N-1)
      SEUIL = 0.001D0
      DF = DFGOAL + SEUIL + 1D0
      CALL IDENTITE(IDENT,N+1)
      CALL DGEMM('T','N',LDQ-2,N+1,LDQ,1D0,Q,LDQ,SIG,LDSIG,
     $              0D0,QSIG,N-1)
      CALL DGEMM('N','N',N-1,N-1,N+1,1D0,QSIG,N-1,Q,LDQ,
     $              0D0,QQ,N-1)
      INDICE = 0
      DO WHILE((ABS(DFGOAL - DF) .gt. SEUIL .and. INDICE .lt. 300))
        PDF = SQRT(A*B)
        DF = 0D0
        CALL ADDITION_MATRICE_CARRE(QQ,N-1,T,LDT,M,N-1,PDF)
        CALL INVERSION(M,N-1,INVERSE,N-1,N-1)
        CALL DGEMM('N','N',LDSIG,N-1,N+1,1D0,SIG,LDSIG,Q,LDQ,
     $              0D0,SIGQ,N+1)
        CALL DGEMM('N','N',N+1,N-1,N-1,1D0,SIGQ,N+1,INVERSE,
     $              N-1,0D0,TEMP,N+1)
        CALL DGEMM('N','T',N+1,N+1,N-1,1D0,TEMP,N+1,Q,
     $              LDQ,0D0,IDENT,N+1)
        DO I = 1,N+1
            DF = DF + (1D0-IDENT(I,I))
        END DO
        IF(DF > DFGOAL) THEN
            B = PDF
        ELSE IF(DF < DFGOAL) THEN 
            A = PDF
        END IF
        INDICE = INDICE + 1
      END DO
      WRITE(*,*) "p -> ",PDF
      RETURN
      END FUNCTION
      

      
