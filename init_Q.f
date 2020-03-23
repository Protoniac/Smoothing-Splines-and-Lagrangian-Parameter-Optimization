      SUBROUTINE INIT_Q(Q,G,LDG)
      IMPLICIT NONE
      INTEGER LDG
      DOUBLE PRECISION G(LDG),Q(LDG+1,LDG-1)
      INTEGER I,J,K
*f2py intent(inplace) q
      K = 0
      DO I=1,LDG-1
        K = K+1
        Q(K,I) = G(K)
        Q(K+1,I) = -1*(G(K)+G(K+1))
        Q(K+2,I) = G(K+1)
        DO J = 1,K-1
            Q(J,I) = 0
        END DO
        DO J = K+3,LDG+1
            Q(J,I) = 0
        END DO
      END DO
      END SUBROUTINE
