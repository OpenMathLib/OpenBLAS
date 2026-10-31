! Coverage driver for the EMPTSW integer-overflow fix in CGESVJ/CGSVJ0
! (issue #1312: EMPTSW = N*(N-1)/2 overflows 32-bit arithmetic for
! N >= 46342, so it is computed in a factored form that never
! overflows).
!
! CGSVJ0 is only reached through CGESVJ's triangular (JOBA='L'/'U')
! blocked path for N > 64, which the xeigtstc SVD tests (JOBA='G')
! never exercise, even with the N=66 pair in TESTING/svd.in.  This
! driver therefore calls CGESVJ with N > 64 lower- and upper-triangular
! matrices (even and odd N), which takes the blocked CGESVJ -> CGSVJ0
! path and covers the fixed lines in both files.  It validates the
! live factorization directly (residual and orthogonality), so there
! is no replica of the EMPTSW expression here to go stale.
program test_cgesvj_emptsw
  use, intrinsic :: iso_fortran_env, only: sp => real32
  implicit none
  integer, parameter :: neven = 66, nodd = 65
  integer, parameter :: lw = 3*neven, lrw = 2*neven
  complex(sp) :: a_even(neven,neven), a0_even(neven,neven)
  complex(sp) :: v_even(neven,neven)
  complex(sp) :: a_odd(nodd,nodd), a0_odd(nodd,nodd)
  complex(sp) :: v_odd(nodd,nodd)
  complex(sp) :: cwork(lw)
  real(sp) :: sva_even(neven), sva_odd(nodd), rwork(lrw)
  real(sp), external :: slamch
  external :: cgesvj
  integer :: i, j, info
  integer :: failures

  failures = 0

  do j = 1, neven
    do i = 1, neven
      if (i == j) then
        a_even(i,j) = cmplx(2.0_sp+0.01_sp*real(i,sp), 0.02_sp*real(i,sp), sp)
      else if (i > j) then
        a_even(i,j) = cmplx(0.1_sp/real(abs(i-j)+1,sp), &
                            0.05_sp/real(abs(i-j)+1,sp), sp)
      else
        a_even(i,j) = cmplx(0.0_sp, 0.0_sp, sp)
      end if
    end do
  end do
  a0_even = a_even
  call cgesvj('L','U','V',neven,neven,a_even,neven,sva_even,neven, &
              v_even,neven,cwork,lw,rwork,lrw,info)
  call check_svd('even', a0_even, a_even, v_even, sva_even, info, &
                 failures)

  do j = 1, nodd
    do i = 1, nodd
      if (i == j) then
        a_odd(i,j) = cmplx(2.0_sp+0.01_sp*real(i,sp), 0.02_sp*real(i,sp), sp)
      else if (i < j) then
        a_odd(i,j) = cmplx(0.1_sp/real(abs(i-j)+1,sp), &
                           0.05_sp/real(abs(i-j)+1,sp), sp)
      else
        a_odd(i,j) = cmplx(0.0_sp, 0.0_sp, sp)
      end if
    end do
  end do
  a0_odd = a_odd
  call cgesvj('U','U','V',nodd,nodd,a_odd,nodd,sva_odd,nodd, &
              v_odd,nodd,cwork,lw,rwork,lrw,info)
  call check_svd('odd ', a0_odd, a_odd, v_odd, sva_odd, info, failures)

  print '(a,i0)', 'CGESVJ EMPTSW failures = ', failures
  if (failures /= 0) stop 1

contains

  subroutine check_svd(tag, a0, u, v, sva, info, failures)
    character(*), intent(in) :: tag
    complex(sp), intent(in) :: a0(:,:), u(:,:), v(:,:)
    real(sp), intent(in) :: sva(:)
    integer, intent(in) :: info
    integer, intent(inout) :: failures
    integer :: n, i, j
    complex(sp), allocatable :: us(:,:), res(:,:), eye(:,:)
    real(sp) :: anorm, rerr, oerr_u, oerr_v, tol, eps

    n = size(sva)
    eps = slamch('Precision')
    tol = real(100*n, sp)*eps
    if (info /= 0) then
      print '(a,a,a,i0)', 'CGESVJ EMPTSW (', tag, ') INFO = ', info
      failures = failures + 1
      return
    end if
    do i = 1, n
      if (sva(i) <= 0.0_sp) failures = failures + 1
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
    eye = (0.0_sp, 0.0_sp)
    do i = 1, n
      eye(i,i) = (1.0_sp, 0.0_sp)
    end do
    oerr_u = maxval(abs(matmul(conjg(transpose(u)), u) - eye))
    oerr_v = maxval(abs(matmul(conjg(transpose(v)), v) - eye))
    print '(a,a,a,3es12.4)', 'CGESVJ EMPTSW (', tag, ') res/orthU/orthV = ', &
          rerr, oerr_u, oerr_v
    if (.not. (rerr <= tol .and. oerr_u <= tol .and. oerr_v <= tol)) &
      failures = failures + 1
  end subroutine

end program
