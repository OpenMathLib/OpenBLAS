x86_64_base_dict = {
    # Level 1 BLAS
    '?rot': {
        's': {
            '_k': 'x86_64/rot_sse.S',
        },
        'd': {
            '_k': 'x86_64/rot_sse2.S',
        },
        'cs': {
            '_k': 'x86_64/zrot_sse.S',
        },
        'zd': {
            '_k': 'x86_64/zrot_sse2.S',
        }
    },
    '?swap': {
        's': {
            '_k': 'x86_64/swap_sse.S',
        },
        'd': {
            '_k': 'x86_64/swap_sse2.S',
        },
        'c': {
            '_k': 'x86_64/zswap_sse.S',
        },
        'z': {
            '_k': 'x86_64/zswap_sse2.S',
        },
    },
    '?scal': {
        's': {
            '_k': 'x86_64/scal_sse.S',
        },
        'd': {
            '_k': 'x86_64/scal_sse2.S',
        },
        'c': {
            '_k': 'x86_64/zscal_sse.S',
        },
        'z': {
            '_k': 'x86_64/zscal_sse2.S',
        },
    },
    '?copy': {
        's': {
            '_k': 'x86_64/copy_sse.S',
        },
        'd': {
            '_k': 'x86_64/copy_sse2.S',
        },
        'c': {
            '_k': 'x86_64/zcopy_sse.S',
        },
        'z': {
            '_k': 'x86_64/zcopy_sse2.S',
        },
    },
    '?axpy': {
        's': {
            '_k': 'x86_64/axpy_sse.S',
        },
        'd': {
            '_k': 'x86_64/axpy_sse2.S',
        },
        'c': {
            '_k': 'x86_64/zaxpy_sse.S',
        },
        'z': {
            '_k': 'x86_64/zaxpy_sse2.S',
        },
    },
    '?axpyc': {
        'c': {
            '_k': 'x86_64/zaxpy_sse.S',
        },
        'z': {
            '_k': 'x86_64/zaxpy_sse2.S',
        },
    },
    '?dot': {
        's': {
            '_k': 'generic/dot.c',
        },
        'd': {
            '_k': 'x86_64/dot_sse2.S',
        },
    },
    '?dotc': {
        'c': {
            '_k': 'x86_64/zdot_sse.S',
        },
        'z': {
            '_k': 'x86_64/zdot_sse2.S',
        },
    },
    '?dotu': {
        'c': {
            '_k': 'x86_64/zdot_sse.S',
        },
        'z': {
            '_k': 'x86_64/zdot_sse2.S',
        },
    },
    '?dsdot': {
        's': {
            '_k': 'generic/dot.c',
        },
        '': {
            '_k': 'generic/dot.c',
        },
    },
    '?nrm2': {
        's': {
            '_k': 'x86_64/nrm2_sse.S',
        },
        'd': {
            '_k': 'x86_64/nrm2.S',
        },
        'c': {
            '_k': 'x86_64/znrm2_sse.S',
        },
        'z': {
            '_k': 'x86_64/znrm2.S',
        },
    },
    '?asum': {
        's': {
            '_k': 'x86_64/asum_sse.S',
        },
        'd': {
            '_k': 'x86_64/asum_sse2.S',
        },
        'c': {
            '_k': 'x86_64/zasum_sse.S',
        },
        'z': {
            '_k': 'x86_64/zasum_sse2.S',
        },
    },
    '?amax': {
        's': {
            '_k': 'x86_64/amax_sse.S',
        },
        'd': {
            '_k': 'x86_64/amax_sse2.S',
        },
        'c': {
            '_k': 'x86_64/zamax_sse.S',
        },
        'z': {
            '_k': 'x86_64/zamax_sse2.S',
        },
    },
    '?sum': {
        's': {
            '_k': 'arm/sum.c',
        },
        'd': {
            '_k': 'arm/sum.c',
        },
        'c': {
            '_k': 'x86_64/zsum_sse.S',
        },
        'z': {
            '_k': 'x86_64/zsum_sse2.S',
        },
    },
    '?amin': {
        's': {
            '_k': 'x86_64/amax_sse.S',
        },
        'd': {
            '_k': 'x86_64/amax_sse2.S',
        },
        'c': {
            '_k': 'x86_64/zamax_sse.S',
        },
        'z': {
            '_k': 'x86_64/zamax_sse2.S',
        },
    },
    'i?amax': {
        's': {
            '_k': 'x86_64/iamax_sse.S',
        },
        'd': {
            '_k': 'x86_64/iamax_sse2.S',
        },
        'c': {
            '_k': 'x86_64/izamax_sse.S',
        },
        'z': {
            '_k': 'x86_64/izamax_sse2.S',
        },
    },
    'i?amin': {
        's': {
            '_k': 'x86_64/iamax_sse.S',
        },
        'd': {
            '_k': 'x86_64/iamax_sse2.S',
        },
        'c': {
            '_k': 'x86_64/izamax_sse.S',
        },
        'z': {
            '_k': 'x86_64/izamax_sse2.S',
        },
    },
    'i?max': {
        's': {
            '_k': 'x86_64/iamax_sse.S',
        },
        'd': {
            '_k': 'x86_64/iamax_sse2.S',
        },
    },
    'i?min': {
        's': {
            '_k': 'x86_64/iamax_sse.S',
        },
        'd': {
            '_k': 'x86_64/iamax_sse2.S',
        },
    },
    '?max': {
        's': {
            '_k': 'x86_64/amax_sse.S',
        },
        'd': {
            '_k': 'x86_64/amax_sse2.S',
        },
    },
    '?min': {
        's': {
            '_k': 'x86_64/amax_sse.S',
        },
        'd': {
            '_k': 'x86_64/amax_sse2.S',
        },
    },
    '?axpby': {
        's': {
            '_k': 'arm/axpby.c',
        },
        'd': {
            '_k': 'arm/axpby.c',
        },
        'c': {
            '_k': 'arm/zaxpby.c',
        },
        'z': {
            '_k': 'arm/zaxpby.c',
        },
    },
    # Level 2 BLAS
    '?gemv': {
        's': {
            '_n': 'x86_64/sgemv_n.c',
            '_t': 'x86_64/sgemv_t.c',
        },
        'd': {
            '_n': 'x86_64/dgemv_n.S',
            '_t': 'x86_64/dgemv_t_4.c',
        },
        'c': {
            '_n': 'x86_64/cgemv_n_4.c',
            '_t': 'x86_64/cgemv_t_4.c',
            '_r': 'x86_64/cgemv_n_4.c',
            '_c': 'x86_64/cgemv_t_4.c',
            '_o': 'x86_64/cgemv_n_4.c',
            '_u': 'x86_64/cgemv_t_4.c',
            '_s': 'x86_64/cgemv_n_4.c',
            '_d': 'x86_64/cgemv_t_4.c',
        },
        'z': {
            '_n': 'x86_64/zgemv_n_4.c',
            '_t': 'x86_64/zgemv_t_4.c',
            '_r': 'x86_64/zgemv_n_4.c',
            '_c': 'x86_64/zgemv_t_4.c',
            '_o': 'x86_64/zgemv_n_4.c',
            '_u': 'x86_64/zgemv_t_4.c',
            '_s': 'x86_64/zgemv_n_4.c',
            '_d': 'x86_64/zgemv_t_4.c',
        },
    },
    '?symv': {
        's': {
            '_U': 'x86_64/symv_U_sse.S',
            '_L': 'x86_64/symv_L_sse.S',
        },
        'd': {
            '_U': 'x86_64/symv_U_sse2.S',
            '_L': 'x86_64/symv_L_sse2.S',
        },
        'c': {
            '_U': 'generic/zsymv_k.c',
            '_L': 'generic/zsymv_k.c',
        },
        'z': {
            '_U': 'x86_64/zsymv_U_sse2.S',
            '_L': 'x86_64/zsymv_L_sse2.S',
        },
    },
    '?lsame': {
        '': {
            '': 'x86_64/lsame.S',
        }
    },
    '?cabs': {
        's': {
            '1': 'x86_64/cabs.S',
        },
        'd': {
            '1': 'x86_64/cabs.S',
        },
    },
    '?gemm3m': {

    },
    '?hemv': {
        'z': {
            '_U': 'x86_64/zsymv_U_sse2.S',
            '_L': 'x86_64/zsymv_L_sse2.S',
        },
    },
    # Level 3 BLAS
    '?gemm_kernel': {
        # done
    },
    '?trmm_kernel': {
        # done
    },
    '?trsm_kernel': {
        # done
    },
    '?gemm': {
        's': {
            '_beta': 'x86_64/gemm_beta.S',
            '_small_matrix_permit': 'generic/gemm_small_matrix_permit.c',
        },
        'd': {
            '_beta': 'x86_64/gemm_beta.S',
        },
        'c': {
            '_beta': 'x86_64/zgemm_beta.S',
        },
        'z': {
            '_beta': 'x86_64/zgemm_beta.S',
        },
    },
    '?trmm': {

    },
    '?hemm': {

    },
    '?trsm': {

    },
    '?symm': {

    },
    '?omatcopy_k': {

    },
    '?imatcopy_k': {

    },
    '?geadd': {

    },
    '?gemm_small_kernel': {

    },
    '?gemm_small_kernel_b0': {

    },
}

subdir('meson_haswell')
subdir('meson_skylakex')
subdir('meson_zen')
subdir('meson_sandybridge')
