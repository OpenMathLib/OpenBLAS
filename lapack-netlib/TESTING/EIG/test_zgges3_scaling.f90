! Check that a negligible nonzero pivot in ZLAQZ0 does not leave fill
! below B's diagonal. Such fill is not rescaled by ZGGES3's final
! ZLASCL('U',...) call and corrupts the decomposition of a small pencil.
program test_zgges3_scaling
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  ! Pad the active 3-by-3 block to enable multishift QZ with the default
  ! IPARMQ crossover NMIN=75, without overriding the tuning routines.
  integer, parameter :: n=80, lwork=200000
  complex(dp) :: a(n,n), b(n,n), a0(n,n), b0(n,n), q(n,n), z(n,n)
  complex(dp) :: alpha(n), beta(n), residual(n,n)
  complex(dp), allocatable :: work(:)
  real(dp) :: rwork(8*n), scale, anorm, bnorm, aerr, berr, lower, tol
  real(dp), external :: dlamch
  logical :: bwork(n), select_none
  external :: select_none, zgges3
  integer :: i, j, info, sdim

  ! Keep the large workspace off the stack, including recursive builds.
  allocate(work(lwork))
  scale=1.0e-280_dp
  a=0; b=0
  do i=1,n
    a(i,i)=1; b(i,i)=1
  end do
  ! A = [0 -1 -1; 1e-20 0 -1; 0 1 -1]
  ! B = [1  1  1; 0 (1+2i)*1e-18 1; 0 0 1]
  a(1,1)=0; a(1,2)=-1; a(1,3)=-1
  a(2,1)=1.0e-20_dp; a(2,2)=0; a(2,3)=-1
  a(3,2)=1; a(3,3)=-1
  b(1,2)=1; b(1,3)=1; b(2,3)=1
  b(2,2)=cmplx(1.0e-18_dp,2.0e-18_dp,dp)
  a=a*scale; b=b*scale
  a0=a; b0=b

  call zgges3('V','V','N',select_none,n,a,n,b,n,sdim,alpha,beta, &
              q,n,z,n,work,lwork,rwork,bwork,info)
  if(info/=0) then
    print '(a,i0)', 'ZGGES3 returned INFO = ',info
    stop 1
  end if

  ! Normalize before reconstruction to avoid underflow in the checks.
  anorm=maxval(abs(a0)); bnorm=maxval(abs(b0))
  residual=a0/anorm-matmul(q,matmul(a/anorm,conjg(transpose(z))))
  aerr=maxval(abs(residual))
  residual=b0/bnorm-matmul(q,matmul(b/bnorm,conjg(transpose(z))))
  berr=maxval(abs(residual))
  lower=0
  do j=1,n
    do i=j+1,n
      lower=max(lower,abs(b(i,j)))
    end do
  end do
  tol=100*n*dlamch('Precision')
  print '(a,es24.16e3)', 'A relative max-entry residual = ',aerr
  print '(a,es24.16e3)', 'B relative max-entry residual = ',berr
  print '(a,es24.16e3)', 'B strict lower triangle maximum = ',lower
  deallocate(work)
  ! Negated comparisons also reject NaNs.
  if(.not.(aerr<=tol .and. berr<=tol .and. lower==0.0_dp)) stop 1
end program

logical function select_none(alpha,beta)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  complex(dp), intent(in) :: alpha,beta
  select_none=.false.
end function
