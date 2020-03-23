        SUBROUTINE CHOLESKY(A,LDA,L,LDL)
        IMPLICIT NONE
        DOUBLE PRECISION A(LDA,*), L(LDL,*), VAL_L
        INTEGER I,J,K,LDA,LDL
*f2py intent(inplace) L
        DO I=1,LDL
            DO J=1,LDL
                IF(J>I-1) THEN
                    VAL_L = 0D0
                    IF (I==J) THEN
                        DO K=1,I-1
                            VAL_L = VAL_L + L(I,K)*L(I,K)
                        END DO
                        L(I,I) = SQRT(A(I,I) - VAL_L)
                    ELSE
                        DO K=1,I-1
                            VAL_L = VAL_L + L(J,K)*L(I,K)
                        END DO
                        L(J,I) = (A(I,J) - VAL_L)/L(I,I)
                    ENDIF
                ELSE 
                    L(J,I) = 0
                ENDIF
            END DO
        END DO
        RETURN
        END SUBROUTINE
