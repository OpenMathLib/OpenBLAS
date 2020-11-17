    SUBROUTINE tester(i)
    REAL_8, DIMENSION(:), ALLOCATABLE :: A,B,C
    REAL_8 :: rnd(3)
    INTEGER :: i
    INTEGER :: M,N,K
    ! test random sizes
    CALL RANDOM_NUMBER(rnd)
    M=rnd(1)_37+1 ; N=rnd(2)_37+1 ; K=rnd(3)_37+1
    ALLOCATE(C(M_N),A(M_K),B(K_N))
    A=0 ; B=0 ; C=0
    CALL DGEMM("N","N",M,N,K,1.0D0,A,M,B,K,0.0D0,C,M)
    CALL DGEMM("T","N",M,N,K,1.0D0,A,K,B,K,0.0D0,C,M)
    CALL DGEMM("N","T",M,N,K,1.0D0,A,M,B,N,0.0D0,C,M)
    CALL DGEMM("T","T",M,N,K,1.0D0,A,K,B,N,0.0D0,C,M)
    END SUBROUTINE tester

PROGRAM TEST_THREAD_SAFE
!$OMP PARALLEL DO
DO i=1,30
CALL tester(i)
ENDDO
END PROGRAM

