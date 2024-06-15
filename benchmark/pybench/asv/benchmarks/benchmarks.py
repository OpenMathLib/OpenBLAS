# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.

'''
class TimeSuite:
    """
    An example benchmark that times the performance of various kinds
    of iterating over dictionaries in Python.
    """
    def setup(self):
        self.d = {}
        for x in range(500):
            self.d[x] = None

    def time_keys(self):
        for key in self.d.keys():
            pass

    def time_values(self):
        for value in self.d.values():
            pass

    def time_range(self):
        d = self.d
        for key in range(500):
            d[key]


class MemSuite:
    def mem_list(self):
        return [0] * 256
'''


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



class Nrm2:

    params = [100, 1000]
    param_names = ["size"]

    def setup(self, n):
        rndm = np.random.RandomState(1234)
        self.x = rndm.uniform(size=(n,)).astype(float)

    def time_dnrm2(self, n):
        run_dnrm2(n, self.x, 1)


# ddot

ddot_sizes = [100, 1000]

def run_ddot(x, y,):
    res = ddot(x, y)
    return res


class DDot:
    params = ddot_sizes
    param_names = ["size"]

    def setup(self, n):
        rndm = np.random.RandomState(1234)
        self.x = np.array(rndm.uniform(size=(n,)), dtype=float)
        self.y = np.array(rndm.uniform(size=(n,)), dtype=float)

    def time_ddot(self, n):
        run_ddot(self.x, self.y)



# daxpy

daxpy_sizes = [100, 1000]

def run_daxpy(x, y,):
    res = daxpy(x, y, a=2.0)
    return res


class Daxpy:
    params = daxpy_sizes
    param_names = ["size"]

    def setup(self, n):
        rndm = np.random.RandomState(1234)
        self.x = np.array(rndm.uniform(size=(n,)), dtype=float)
        self.y = np.array(rndm.uniform(size=(n,)), dtype=float)

    def time_daxpy(self, n):
        run_daxpy(self.x, self.y)



# ### BLAS level 3 ###

# dgemm

gemm_sizes = [100, 1000]

def run_dgemm(a, b, c):
    alpha = 1.0
    res = dgemm(alpha, a, b, c=c, overwrite_c=True)
    return res


class Dgemm:
    params = gemm_sizes
    param_names = ["size"]

    def setup(self, n):
        rndm = np.random.RandomState(1234)
        self.a = np.array(rndm.uniform(size=(n, n)), dtype=float, order='F')
        self.b = np.array(rndm.uniform(size=(n, n)), dtype=float, order='F')
        self.c = np.empty((n, n), dtype=float, order='F')

    def time_dgemm(self, n):
        run_dgemm(self.a, self.b, self.c)


# dsyrk

syrk_sizes = [100, 1000]


def run_dsyrk(a, c):
    res = dsyrk(1.0, a, c=c, overwrite_c=True)
    return res


class DSyrk:
    params = syrk_sizes
    param_names = ["size"]

    def setup(self, n):
        rndm = np.random.RandomState(1234)
        self.a = np.array(rndm.uniform(size=(n, n)), dtype=float, order='F')
        self.c = np.empty((n, n), dtype=float, order='F')

    def time_dsyrk(self, n):
        run_dsyrk(self.a, self.c)


# ### LAPACK ###

# linalg.solve

dgesv_sizes = [100, 1000]


def run_dgesv(a, b):
    res = dgesv(a, b, overwrite_a=True, overwrite_b=True)
    return res


class Dgesv:
    params = dgesv_sizes
    param_names = ["size"]

    def setup(self, n):
        rndm = np.random.RandomState(1234)
        self.a = (np.array(rndm.uniform(size=(n, n)), dtype=float, order='F') +
                  np.eye(n, order='F'))
        self.b = np.array(rndm.uniform(size=(n, 1)), order='F')

    def time_dgesv(self, n):
        run_dgesv(self.a, self.b)

      # XXX: how to run asserts?
      #  lu, piv, x, info = benchmark(run_gesv, a, b)
      #  assert lu is a
      #  assert x is b
      #  assert info == 0


# linalg.svd

dgesdd_sizes = ["100, 5", "1000, 222"]


def run_dgesdd(a, lwork):
    res = dgesdd(a, lwork=lwork, full_matrices=False, overwrite_a=False)
    return res


class Dgesdd:
    params = dgesdd_sizes
    param_names = ["(m, n)"]

    def setup(self, mn):
        m, n = (int(x) for x in mn.split(","))

        rndm = np.random.RandomState(1234)
        a = np.array(rndm.uniform(size=(m, n)), dtype=float, order='F')

        lwork, info = dgesdd_lwork(m, n)
        lwork = int(lwork)
        assert info == 0

        self.a, self.lwork = a, lwork

    def time_dgesdd(self, mn):
        run_dgesdd(self.a, self.lwork)


# linalg.eigh

dsyev_sizes = [50, 200]


def run_dsyev(a, lwork):
    res = dsyev(a, lwork=lwork, overwrite_a=True)
    return res


class Dsyev:
    params = dsyev_sizes
    param_names = ["size"]

    def setup(self, n):
        rndm = np.random.RandomState(1234)
        a = rndm.uniform(size=(n, n))
        a = np.asarray(a + a.T, dtype=float, order='F')
        a_ = a.copy()

        lwork, info = dsyev_lwork(n)
        lwork = int(lwork)
        assert info == 0

        self.a = a_
        self.lwork = lwork

    def time_dsyev(self, n):
        run_dsyev(self.a, self.lwork)

