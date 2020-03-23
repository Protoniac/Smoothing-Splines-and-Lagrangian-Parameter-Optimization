      SUBROUTINE INIT_T(T,H,LDH)
      IMPLICIT NONE
      INTEGER LDH
      DOUBLE PRECISION T(LDH-1,*)
      DOUBLE PRECISION H(LDH)
      INTEGER I,J
*f2py intent(inplace) t
      DO I=1,(LDH-1)
        DO J=1,(LDH-1)
            IF(I==J) THEN
                T(I,J) = 2*(H(I)+H(I+1))/3
            ELSE IF(I+1 == J) THEN 
                T(I,J) = H(J)/3
            ELSE IF(I-1 == J) THEN
                T(I,J) = H(I)/3
            ELSE
                T(I,J) = 0
            END IF
        END DO
      END DO
      END SUBROUTINE
