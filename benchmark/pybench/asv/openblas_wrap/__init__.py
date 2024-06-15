"""
Trampoline to hide the LAPACK details (scipy.lapack.linalg or scipy_openblas32 or...)
from benchmarking.
"""

__version__ = "0.1"  

import scipy_openblas32   # preload symbols. typically done in _distributor_init.py

#from scipy.linalg.blas import (
from ._flapack import (
    # level 1
    scipy_dnrm2 as dnrm2,
    scipy_ddot as ddot,
    scipy_daxpy as daxpy,
    # level 3
    scipy_dgemm as dgemm,
    scipy_dsyrk as dsyrk,
)

#from scipy.linalg.lapack import (
from openblas_wrap._flapack import (
    # linalg.solve
    scipy_dgesv as dgesv,
    # linalg.svd
    scipy_dgesdd as dgesdd, scipy_dgesdd_lwork as dgesdd_lwork,
    # linalg.eigh
    scipy_dsyev as dsyev, scipy_dsyev_lwork as dsyev_lwork
)
