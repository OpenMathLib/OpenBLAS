"""
Trampoline to hide the LAPACK details (scipy.lapack.linalg or scipy_openblas32 or...)
from benchmarking.
"""

__version__ = "0.1"  


#from scipy.linalg.blas import (
from ._flapack import (
    # level 1
    dnrm2 as dnrm2,
    ddot as ddot,
    daxpy as daxpy,
    # level 3
    dgemm as dgemm,
    dsyrk as dsyrk,
)

#from scipy.linalg.lapack import (
from openblas_wrap._flapack import (
    # linalg.solve
    dgesv as dgesv,
    # linalg.svd
    dgesdd as dgesdd, dgesdd_lwork as dgesdd_lwork,
    # linalg.eigh
    dsyev as dsyev, dsyev_lwork as dsyev_lwork
)
