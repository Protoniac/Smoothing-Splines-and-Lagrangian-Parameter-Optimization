      SUBROUTINE INIT_HG(H,G,X,LDX)
      IMPLICIT NONE
      INTEGER LDX
      DOUBLE PRECISION H(LDX-1),G(LDX-1),X(LDX)
      INTEGER I
*f2py intent(inplace) h
*f2py intent(inplace) g
      DO I=1,(LDX-1)
        H(I) = X(I+1) - X(I)
        G(I) = 1/H(I)
      END DO
      END SUBROUTINE
