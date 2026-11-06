! Check that a negligible nonzero pivot in SLAQZ0 does not leave fill
! below B's diagonal before SGGES3 restores the original scale.
program test_sgges3_scaling
  use, intrinsic :: iso_fortran_env, only: sp => real32
  implicit none
  integer, parameter :: n=80, lwork=200000
  real(sp) :: a(n,n), b(n,n), a0(n,n), b0(n,n), q(n,n), z(n,n)
  real(sp) :: alphar(n), alphai(n), beta(n), residual(n,n)
  real(sp), allocatable :: work(:)
  real(sp) :: scale, anorm, bnorm, aerr, berr, lower, tol
  real(sp), external :: slamch
  logical :: bwork(n), select_none
  external :: select_none, sgges3
  integer :: i, j, info, sdim

  allocate(work(lwork))
  scale=1.0e-30_sp
  a=0; b=0
  do i=1,n
    a(i,i)=1; b(i,i)=1
  end do
  a(1,1)=0; a(1,2)=-1; a(1,3)=-1
  a(2,1)=1.0e-10_sp; a(2,2)=0; a(2,3)=-1
  a(3,2)=1; a(3,3)=-1
  b(1,2)=1; b(1,3)=1; b(2,3)=1
  b(2,2)=1.0e-8_sp
  a=a*scale; b=b*scale
  a0=a; b0=b

  call sgges3('V','V','N',select_none,n,a,n,b,n,sdim,alphar,alphai, &
              beta,q,n,z,n,work,lwork,bwork,info)
  if(info/=0) then
    print '(a,i0)', 'SGGES3 returned INFO = ',info
    stop 1
  end if
  anorm=maxval(abs(a0)); bnorm=maxval(abs(b0))
  residual=a0/anorm-matmul(q,matmul(a/anorm,transpose(z)))
  aerr=maxval(abs(residual))
  residual=b0/bnorm-matmul(q,matmul(b/bnorm,transpose(z)))
  berr=maxval(abs(residual))
  lower=0
  do j=1,n
    do i=j+1,n
      lower=max(lower,abs(b(i,j)))
    end do
  end do
  tol=100*n*slamch('Precision')
  print '(a,es24.16e3)', 'A relative max-entry residual = ',aerr
  print '(a,es24.16e3)', 'B relative max-entry residual = ',berr
  print '(a,es24.16e3)', 'B strict lower triangle maximum = ',lower
  deallocate(work)
  if(.not.(aerr<=tol .and. berr<=tol .and. lower==0.0_sp)) stop 1
end program

logical function select_none(alphar,alphai,beta)
  use, intrinsic :: iso_fortran_env, only: sp => real32
  implicit none
  real(sp), intent(in) :: alphar,alphai,beta
  select_none=.false.
end function
