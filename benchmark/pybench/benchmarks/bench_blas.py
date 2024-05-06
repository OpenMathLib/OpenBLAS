import pytest
import numpy as np
from openblas_wrap import (
    # level 1
    dnrm2, ddot, daxpy,
    # level 3
    dgemm, dsyrk,
    # lapack
    dgesv,                   # linalg.solve
    dgesdd, dgesdd_lwork,    # linalg.svd
    dsyev, dsyev_lwork,      # linalg.eigh
)

# ### BLAS level 1 ###

# dnrm2

dnrm2_sizes = [100, 1000]

def run_dnrm2(n, x, incx):
    res = dnrm2(x, n, incx=incx)
    return res


@pytest.mark.parametrize('n', dnrm2_sizes)
def test_nrm2(benchmark, n):
    rndm = np.random.RandomState(1234)
    x = np.array(rndm.uniform(size=(n,)), dtype=float)
    result = benchmark(run_dnrm2, n, x, 1)


# ddot

ddot_sizes = [100, 1000]

def run_ddot(x, y,):
    res = ddot(x, y)
    return res


@pytest.mark.parametrize('n', ddot_sizes)
def test_dot(benchmark, n):
    rndm = np.random.RandomState(1234)
    x = np.array(rndm.uniform(size=(n,)), dtype=float)
    y = np.array(rndm.uniform(size=(n,)), dtype=float)
    result = benchmark(run_ddot, x, y)


# daxpy

daxpy_sizes = [100, 1000]

def run_daxpy(x, y,):
    res = daxpy(x, y, a=2.0)
    return res


@pytest.mark.parametrize('n', daxpy_sizes)
def test_daxpy(benchmark, n):
    rndm = np.random.RandomState(1234)
    x = np.array(rndm.uniform(size=(n,)), dtype=float)
    y = np.array(rndm.uniform(size=(n,)), dtype=float)
    result = benchmark(run_daxpy, x, y)




# ### BLAS level 3 ###

# dgemm

gemm_sizes = [100, 1000]

def run_gemm(a, b, c):
    alpha = 1.0
    res = dgemm(alpha, a, b, c=c, overwrite_c=True)
    return res


@pytest.mark.parametrize('n', gemm_sizes)
def test_gemm(benchmark, n):
    rndm = np.random.RandomState(1234)
    a = np.array(rndm.uniform(size=(n, n)), dtype=float, order='F')
    b = np.array(rndm.uniform(size=(n, n)), dtype=float, order='F')
    c = np.empty((n, n), dtype=float, order='F')
    result = benchmark(run_gemm, a, b, c)
    assert result is c


# dsyrk

syrk_sizes = [100, 1000]


def run_syrk(a, c):
    res = dsyrk(1.0, a, c=c, overwrite_c=True)
    return res


@pytest.mark.parametrize('n', syrk_sizes)
def test_syrk(benchmark, n):
    rndm = np.random.RandomState(1234)
    a = np.array(rndm.uniform(size=(n, n)), dtype=float, order='F')
    c = np.empty((n, n), dtype=float, order='F')
    result = benchmark(run_syrk, a, c)
    assert result is c


# ### LAPACK ###

# linalg.solve

gesv_sizes = [100, 1000]


def run_gesv(a, b):
    res = dgesv(a, b, overwrite_a=True, overwrite_b=True)
    return res


@pytest.mark.parametrize('n', gesv_sizes)
def test_gesv(benchmark, n):
    rndm = np.random.RandomState(1234)
    a = (np.array(rndm.uniform(size=(n, n)), dtype=float, order='F') +
         np.eye(n, order='F'))
    b = np.array(rndm.uniform(size=(n, 1)), order='F')
    lu, piv, x, info = benchmark(run_gesv, a, b)
    assert lu is a
    assert x is b
    assert info == 0


# linalg.svd

gesdd_sizes = [(100, 5), (1000, 222)]


def run_gesdd(a, lwork):
    res = dgesdd(a, lwork=lwork, full_matrices=False, overwrite_a=False)
    return res


@pytest.mark.parametrize('mn', gesdd_sizes)
def test_gesdd(benchmark, mn):
    m, n = mn
    rndm = np.random.RandomState(1234)
    a = np.array(rndm.uniform(size=(m, n)), dtype=float, order='F')

    lwork, info = dgesdd_lwork(m, n)
    lwork = int(lwork)
    assert info == 0

    u, s, vt, info = benchmark(run_gesdd, a, lwork)

    assert info == 0
    np.testing.assert_allclose(u @ np.diag(s) @ vt, a, atol=1e-13)


# linalg.eigh

syev_sizes = [50, 200]


def run_syev(a, lwork):
    res = dsyev(a, lwork=lwork, overwrite_a=True)
    return res


@pytest.mark.parametrize('n', syev_sizes)
def test_syev(benchmark, n):
    rndm = np.random.RandomState(1234)
    a = rndm.uniform(size=(n, n))
    a = np.asarray(a + a.T, dtype=float, order='F')
    a_ = a.copy()

    lwork, info = dsyev_lwork(n)
    lwork = int(lwork)
    assert info == 0

    w, v, info = benchmark(run_syev, a, lwork)

    assert info == 0
    assert a is v  # overwrite_a=True


