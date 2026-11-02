! Coverage driver for the EMPTSW integer-overflow fix in SGESVJ/SGSVJ0
! (issue #1312: EMPTSW = N*(N-1)/2 overflows 32-bit arithmetic for
! N >= 46342, so it is computed in a factored form that never
! overflows).
!
! SGSVJ0 is only reached through SGESVJ's triangular (JOBA='L'/'U')
! blocked path for N > 64, which the xeigtsts SVD tests (JOBA='G')
! never exercise, even with the N=66 pair in TESTING/svd.in.  This
! driver therefore calls SGESVJ with N > 64 lower- and upper-triangular
! matrices (even and odd N), which takes the blocked SGESVJ -> SGSVJ0
! path and covers the fixed lines in both files.  It validates the
! live factorization directly (residual and orthogonality), so there
! is no replica of the EMPTSW expression here to go stale.
program test_sgesvj_emptsw
  use, intrinsic :: iso_fortran_env, only: sp => real32
  implicit none
  integer, parameter :: neven = 66, nodd = 65, lw = 3*neven
  real(sp) :: a_even(neven,neven), a0_even(neven,neven)
  real(sp) :: v_even(neven,neven), sva_even(neven)
  real(sp) :: a_odd(nodd,nodd), a0_odd(nodd,nodd)
  real(sp) :: v_odd(nodd,nodd), sva_odd(nodd)
  real(sp) :: work(lw)
  real(sp), external :: slamch
  external :: sgesvj
  integer :: i, j, info
  integer :: failures

  failures = 0

  do j = 1, neven
    do i = 1, neven
      if (i == j) then
        a_even(i,j) = 2.0_sp + 0.01_sp*real(i,sp)
      else if (i > j) then
        a_even(i,j) = 0.1_sp/real(abs(i-j)+1,sp)
      else
        a_even(i,j) = 0.0_sp
      end if
    end do
  end do
  a0_even = a_even
  call sgesvj('L','U','V',neven,neven,a_even,neven,sva_even,neven, &
              v_even,neven,work,lw,info)
  call check_svd('even', a0_even, a_even, v_even, sva_even, info, &
                 failures)

  do j = 1, nodd
    do i = 1, nodd
      if (i == j) then
        a_odd(i,j) = 2.0_sp + 0.01_sp*real(i,sp)
      else if (i < j) then
        a_odd(i,j) = 0.1_sp/real(abs(i-j)+1,sp)
      else
        a_odd(i,j) = 0.0_sp
      end if
    end do
  end do
  a0_odd = a_odd
  call sgesvj('U','U','V',nodd,nodd,a_odd,nodd,sva_odd,nodd, &
              v_odd,nodd,work,lw,info)
  call check_svd('odd ', a0_odd, a_odd, v_odd, sva_odd, info, failures)

  print '(a,i0)', 'SGESVJ EMPTSW failures = ', failures
  if (failures /= 0) stop 1

contains

  subroutine check_svd(tag, a0, u, v, sva, info, failures)
    character(*), intent(in) :: tag
    real(sp), intent(in) :: a0(:,:), u(:,:), v(:,:)
    real(sp), intent(in) :: sva(:)
    integer, intent(in) :: info
    integer, intent(inout) :: failures
    integer :: n, i, j
    real(sp), allocatable :: us(:,:), res(:,:), eye(:,:)
    real(sp) :: anorm, rerr, oerr_u, oerr_v, tol, eps

    n = size(sva)
    eps = slamch('Precision')
    tol = real(100*n, sp)*eps
    if (info /= 0) then
      print '(a,a,a,i0)', 'SGESVJ EMPTSW (', tag, ') INFO = ', info
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
    res = a0 - matmul(us, transpose(v))
    rerr = maxval(abs(res)) / anorm
    eye = 0.0_sp
    do i = 1, n
      eye(i,i) = 1.0_sp
    end do
    oerr_u = maxval(abs(matmul(transpose(u), u) - eye))
    oerr_v = maxval(abs(matmul(transpose(v), v) - eye))
    print '(a,a,a,3es12.4)', 'SGESVJ EMPTSW (', tag, ') res/orthU/orthV = ', &
          rerr, oerr_u, oerr_v
    if (.not. (rerr <= tol .and. oerr_u <= tol .and. oerr_v <= tol)) &
      failures = failures + 1
  end subroutine

end program
