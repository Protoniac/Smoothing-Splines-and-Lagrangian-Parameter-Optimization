!    -*- f90 -*-
! Note: the context of this file is case sensitive.

subroutine cholesky(a,lda,l,ldl) ! in cholesky.f
    double precision dimension(lda,*) :: a
    integer, optional,check(shape(a,0)==lda),depend(a) :: lda=shape(a,0)
    double precision dimension(ldl,*),intent(inplace) :: l
    integer, optional,check(shape(l,0)==ldl),depend(l) :: ldl=shape(l,0)
end subroutine cholesky
subroutine init_t(t,h,ldh) ! in init_T.f
    double precision dimension(ldh - 1,*),intent(inplace) :: t
    double precision dimension(ldh),depend(ldh) :: h
    integer, optional,check((shape(t,0)+1)==ldh),depend(t) :: ldh=(shape(t,0)+1)
end subroutine init_t
subroutine init_hg(h,g,x,ldx) ! in init_HG.f
    double precision dimension(ldx - 1),intent(inplace) :: h
    double precision dimension(ldx - 1),intent(inplace),depend(ldx) :: g
    double precision dimension(ldx),depend(ldx) :: x
    integer, optional,check((len(h)+1)>=ldx),depend(h) :: ldx=(len(h)+1)
end subroutine init_hg
subroutine init_q(q,g,ldg) ! in init_Q.f
    double precision dimension(ldg + 1,ldg - 1),intent(inplace) :: q
    double precision dimension(ldg),depend(ldg) :: g
    integer, optional,check((shape(q,0)-1)==ldg),depend(q) :: ldg=(shape(q,0)-1)
end subroutine init_q
subroutine calcul_a(q,ldq,c,ldc,y,ldy,a,lda,p) ! in outils.f
    double precision dimension(ldq,*) :: q
    integer, optional,check(shape(q,0)==ldq),depend(q) :: ldq=shape(q,0)
    double precision dimension(ldc) :: c
    integer, optional,check(len(c)>=ldc),depend(c) :: ldc=len(c)
    double precision dimension(ldy) :: y
    integer, optional,check(len(y)>=ldy),depend(y) :: ldy=len(y)
    double precision dimension(lda) :: a
    integer, optional,check(len(a)>=lda),depend(a) :: lda=len(a)
    double precision :: p
end subroutine calcul_a
subroutine calcul_d(c,ldc,h,ldh,d,ldd) ! in outils.f
    double precision dimension(ldc) :: c
    integer, optional,check(len(c)>=ldc),depend(c) :: ldc=len(c)
    double precision dimension(ldh) :: h
    integer, optional,check(len(h)>=ldh),depend(h) :: ldh=len(h)
    double precision dimension(ldd) :: d
    integer, optional,check(len(d)>=ldd),depend(d) :: ldd=len(d)
end subroutine calcul_d
subroutine calcul_b(a,lda,c,ldc,d,ldd,h,ldh,b,ldb) ! in outils.f
    double precision dimension(lda) :: a
    integer, optional,check(len(a)>=lda),depend(a) :: lda=len(a)
    double precision dimension(ldc) :: c
    integer, optional,check(len(c)>=ldc),depend(c) :: ldc=len(c)
    double precision dimension(ldd) :: d
    integer, optional,check(len(d)>=ldd),depend(d) :: ldd=len(d)
    double precision dimension(ldh) :: h
    integer, optional,check(len(h)>=ldh),depend(h) :: ldh=len(h)
    double precision dimension(ldb) :: b
    integer, optional,check(len(b)>=ldb),depend(b) :: ldb=len(b)
end subroutine calcul_b
subroutine addition_matrice_carre(a,lda,b,ldb,m,ldm,alpha) ! in outils.f
    double precision dimension(lda,*) :: a
    integer, optional,check(shape(a,0)==lda),depend(a) :: lda=shape(a,0)
    double precision dimension(ldb,*) :: b
    integer, optional,check(shape(b,0)==ldb),depend(b) :: ldb=shape(b,0)
    double precision dimension(ldm,*) :: m
    integer, optional,check(shape(m,0)==ldm),depend(m) :: ldm=shape(m,0)
    double precision :: alpha
end subroutine addition_matrice_carre
subroutine descente(u,ldu,b,ldb) ! in outils.f
    double precision dimension(ldu,*) :: u
    integer, optional,check(shape(u,0)==ldu),depend(u) :: ldu=shape(u,0)
    double precision dimension(ldb) :: b
    integer, optional,check(len(b)>=ldb),depend(b) :: ldb=len(b)
end subroutine descente
subroutine remontee(u,ldu,b,ldb) ! in outils.f
    double precision dimension(ldu,*) :: u
    integer, optional,check(shape(u,0)==ldu),depend(u) :: ldu=shape(u,0)
    double precision dimension(ldb) :: b
    integer, optional,check(len(b)>=ldb),depend(b) :: ldb=len(b)
end subroutine remontee
subroutine identite(a,lda) ! in outils.f
    double precision dimension(lda,lda) :: a
    integer, optional,check(shape(a,0)==lda),depend(a) :: lda=shape(a,0)
end subroutine identite
subroutine inversion(a,lda,inverse,ldinverse,ldcinverse) ! in outils.f
    double precision dimension(lda,*) :: a
    integer, optional,check(shape(a,0)==lda),depend(a) :: lda=shape(a,0)
    double precision dimension(ldinverse,ldcinverse) :: inverse
    integer, optional,check(shape(inverse,0)==ldinverse),depend(inverse) :: ldinverse=shape(inverse,0)
    integer, optional,check(shape(inverse,1)==ldcinverse),depend(inverse) :: ldcinverse=shape(inverse,1)
end subroutine inversion
function norme2(u,ldu,v,ldv) ! in outils.f
    double precision dimension(ldu) :: u
    integer, optional,check(len(u)>=ldu),depend(u) :: ldu=len(u)
    double precision dimension(ldv) :: v
    integer, optional,check(len(v)>=ldv),depend(v) :: ldv=len(v)
    double precision :: norme2
end function norme2
function newton_p(q,ldq,t,ldt,sig,ldsig,y,ldy,n,s,iteration) ! in newton.f
    double precision dimension(ldq,*) :: q
    integer, optional,check(shape(q,0)==ldq),depend(q) :: ldq=shape(q,0)
    double precision dimension(ldt,*) :: t
    integer, optional,check(shape(t,0)==ldt),depend(t) :: ldt=shape(t,0)
    double precision dimension(ldsig,*) :: sig
    integer, optional,check(shape(sig,0)==ldsig),depend(sig) :: ldsig=shape(sig,0)
    double precision dimension(ldy) :: y
    integer, optional,check(len(y)>=ldy),depend(y) :: ldy=len(y)
    integer :: n
    double precision :: s
    integer :: iteration
    double precision :: newton_p
end function newton_p
subroutine dfp(q,ldq,t,ldt,sig,ldsig,n,res,ldres,ind,ldind,a,b,h) ! in df.f
    double precision dimension(ldq,*) :: q
    integer, optional,check(shape(q,0)==ldq),depend(q) :: ldq=shape(q,0)
    double precision dimension(ldt,*) :: t
    integer, optional,check(shape(t,0)==ldt),depend(t) :: ldt=shape(t,0)
    double precision dimension(ldsig,*) :: sig
    integer, optional,check(shape(sig,0)==ldsig),depend(sig) :: ldsig=shape(sig,0)
    integer :: n
    double precision dimension(ldres),intent(inplace) :: res
    integer, optional,check(len(res)>=ldres),depend(res) :: ldres=len(res)
    double precision dimension(ldind),intent(inplace) :: ind
    integer, optional,check(len(ind)>=ldind),depend(ind) :: ldind=len(ind)
    double precision :: a
    double precision :: b
    double precision :: h
end subroutine dfp
function pdf(q,ldq,t,ldt,sig,ldsig,n,dfgoal,a,b) ! in df.f
    double precision dimension(ldq,*) :: q
    integer, optional,check(shape(q,0)==ldq),depend(q) :: ldq=shape(q,0)
    double precision dimension(ldt,*) :: t
    integer, optional,check(shape(t,0)==ldt),depend(t) :: ldt=shape(t,0)
    double precision dimension(ldsig,*) :: sig
    integer, optional,check(shape(sig,0)==ldsig),depend(sig) :: ldsig=shape(sig,0)
    integer :: n
    double precision :: dfgoal
    double precision :: a
    double precision :: b
    double precision :: pdf
end function pdf
subroutine splines(y,ldy,h,ldh,q,ldq,t,ldt,sig,ldsig,n,p,a,lda,b,ldb,c,ldc,d,ldd) ! in splines.f
    double precision dimension(ldy) :: y
    integer, optional,check(len(y)>=ldy),depend(y) :: ldy=len(y)
    double precision dimension(ldh) :: h
    integer, optional,check(len(h)>=ldh),depend(h) :: ldh=len(h)
    double precision dimension(ldq,*) :: q
    integer, optional,check(shape(q,0)==ldq),depend(q) :: ldq=shape(q,0)
    double precision dimension(ldt,*) :: t
    integer, optional,check(shape(t,0)==ldt),depend(t) :: ldt=shape(t,0)
    double precision dimension(ldsig,*) :: sig
    integer, optional,check(shape(sig,0)==ldsig),depend(sig) :: ldsig=shape(sig,0)
    integer :: n
    double precision :: p
    double precision dimension(lda),intent(inplace) :: a
    integer, optional,check(len(a)>=lda),depend(a) :: lda=len(a)
    double precision dimension(ldb),intent(inplace) :: b
    integer, optional,check(len(b)>=ldb),depend(b) :: ldb=len(b)
    double precision dimension(ldc),intent(inplace) :: c
    integer, optional,check(len(c)>=ldc),depend(c) :: ldc=len(c)
    double precision dimension(ldd),intent(inplace) :: d
    integer, optional,check(len(d)>=ldd),depend(d) :: ldd=len(d)
end subroutine splines

! This file was auto-generated with f2py (version:2).
! See http://cens.ioc.ee/projects/f2py2e/
