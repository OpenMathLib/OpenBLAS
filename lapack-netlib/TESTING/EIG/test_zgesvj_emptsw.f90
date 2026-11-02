! Coverage driver for the EMPTSW integer-overflow fix in ZGESVJ/ZGSVJ0
! (issue #1312: EMPTSW = N*(N-1)/2 overflows 32-bit arithmetic for
! N >= 46342, so it is computed in a factored form that never
! overflows).
!
! ZGSVJ0 is only reached through ZGESVJ's triangular (JOBA='L'/'U')
! blocked path for N > 64, which the xeigtstz SVD tests (JOBA='G')
! never exercise, even with the N=66 pair in TESTING/svd.in.  This
! driver therefore calls ZGESVJ with N > 64 lower- and upper-triangular
! matrices (even and odd N), which takes the blocked ZGESVJ -> ZGSVJ0
! path and covers the fixed lines in both files.  It validates the
! live factorization directly (residual and orthogonality), so there
! is no replica of the EMPTSW expression here to go stale.
program test_zgesvj_emptsw
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  integer, parameter :: neven = 66, nodd = 65
  integer, parameter :: lw = 3*neven, lrw = 2*neven
  complex(dp) :: a_even(neven,neven), a0_even(neven,neven)
  complex(dp) :: v_even(neven,neven)
  complex(dp) :: a_odd(nodd,nodd), a0_odd(nodd,nodd)
  complex(dp) :: v_odd(nodd,nodd)
  complex(dp) :: cwork(lw)
  real(dp) :: sva_even(neven), sva_odd(nodd), rwork(lrw)
  real(dp), external :: dlamch
  external :: zgesvj
  integer :: i, j, info
  integer :: failures

  failures = 0

  do j = 1, neven
    do i = 1, neven
      if (i == j) then
        a_even(i,j) = cmplx(2.0_dp+0.01_dp*real(i,dp), 0.02_dp*real(i,dp), dp)
      else if (i > j) then
        a_even(i,j) = cmplx(0.1_dp/real(abs(i-j)+1,dp), &
                            0.05_dp/real(abs(i-j)+1,dp), dp)
      else
        a_even(i,j) = cmplx(0.0_dp, 0.0_dp, dp)
      end if
    end do
  end do
  a0_even = a_even
  call zgesvj('L','U','V',neven,neven,a_even,neven,sva_even,neven, &
              v_even,neven,cwork,lw,rwork,lrw,info)
  call check_svd('even', a0_even, a_even, v_even, sva_even, info, &
                 failures)

  do j = 1, nodd
    do i = 1, nodd
      if (i == j) then
        a_odd(i,j) = cmplx(2.0_dp+0.01_dp*real(i,dp), 0.02_dp*real(i,dp), dp)
      else if (i < j) then
        a_odd(i,j) = cmplx(0.1_dp/real(abs(i-j)+1,dp), &
                           0.05_dp/real(abs(i-j)+1,dp), dp)
      else
        a_odd(i,j) = cmplx(0.0_dp, 0.0_dp, dp)
      end if
    end do
  end do
  a0_odd = a_odd
  call zgesvj('U','U','V',nodd,nodd,a_odd,nodd,sva_odd,nodd, &
              v_odd,nodd,cwork,lw,rwork,lrw,info)
  call check_svd('odd ', a0_odd, a_odd, v_odd, sva_odd, info, failures)

  print '(a,i0)', 'ZGESVJ EMPTSW failures = ', failures
  if (failures /= 0) stop 1

contains

  subroutine check_svd(tag, a0, u, v, sva, info, failures)
    character(*), intent(in) :: tag
    complex(dp), intent(in) :: a0(:,:), u(:,:), v(:,:)
    real(dp), intent(in) :: sva(:)
    integer, intent(in) :: info
    integer, intent(inout) :: failures
    integer :: n, i, j
    complex(dp), allocatable :: us(:,:), res(:,:), eye(:,:)
    real(dp) :: anorm, rerr, oerr_u, oerr_v, tol, eps

    n = size(sva)
    eps = dlamch('Precision')
    tol = real(100*n, dp)*eps
    if (info /= 0) then
      print '(a,a,a,i0)', 'ZGESVJ EMPTSW (', tag, ') INFO = ', info
      failures = failures + 1
      return
    end if
    do i = 1, n
      if (sva(i) <= 0.0_dp) failures = failures + 1
      if (i > 1) then
        if (sva(i-1) < sva(i)) failures = failures + 1
      end if
    end do
    allocate(us(n,n), res(n,n), eye(n,n))
    anorm = maxval(abs(a0))
    do j = 1, n
      do i = 1, n
        us(i,j) = u(i,j)*sva(j)
      end do
    end do
    res = a0 - matmul(us, conjg(transpose(v)))
    rerr = maxval(abs(res)) / anorm
    eye = (0.0_dp, 0.0_dp)
    do i = 1, n
      eye(i,i) = (1.0_dp, 0.0_dp)
    end do
    oerr_u = maxval(abs(matmul(conjg(transpose(u)), u) - eye))
    oerr_v = maxval(abs(matmul(conjg(transpose(v)), v) - eye))
    print '(a,a,a,3es12.4)', 'ZGESVJ EMPTSW (', tag, ') res/orthU/orthV = ', &
          rerr, oerr_u, oerr_v
    if (.not. (rerr <= tol .and. oerr_u <= tol .and. oerr_v <= tol)) &
      failures = failures + 1
  end subroutine

end program
