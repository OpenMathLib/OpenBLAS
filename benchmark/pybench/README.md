# Continuous benchmarking of OpenBLAS performance

We run a set of benchmarks of subset of OpenBLAS functionality.

## Benchmark runner

[![CodSpeed Badge](https://img.shields.io/endpoint?url=https://codspeed.io/badge.json)](https://codspeed.io/OpenMathLib/OpenBLAS/)

Click on [benchmarks](https://codspeed.io/OpenMathLib/OpenBLAS/benchmarks) to see the performance of a particular benchmark over time;
Click on [branches](https://codspeed.io/OpenMathLib/OpenBLAS/branches/) and then on the last PR link to see the flamegraphs.

## What are the benchmarks

We run raw BLAS/LAPACK subroutines, via f2py-generated python wrappers. The wrappers themselves are equivalent to [those from SciPy](https://docs.scipy.org/doc/scipy/reference/linalg.lapack.html).
In fact, the wrappers _are_ from SciPy, we take a small subset simply to avoid having to build the whole SciPy for each CI run.


## Adding a new benchmark

`.github/workflows/codspeed-bench.yml` does all the orchestration on CI.

Benchmarks live in the `benchmark/pybench` directory. It is organized as follows:

- benchmarks themselves live in the `benchmarks` folder. Note that the LAPACK routines are imported from the `openblas_wrap` package.
- the `openblas_wrap` package is a simple trampoline: it contains an f2py extension, `_flapack`, which talks to OpenBLAS, and exports the python names in its `__init__.py`.
This way, the `openblas_wrap` package shields the benchmarks from the details of where a particular LAPACK function comes from. If wanted, you may for instance swap the `_flapack` extension to
`scipy.linalg.blas` and `scipy.linalg.lapack`.

To change parameters of an existing benchmark, edit python files in the `benchmark/pybench/benchmarks` directory.

To add a benchmark for a new BLAS or LAPACK function, you need to:

- add an f2py wrapper for the bare LAPACK function. You can simply copy a wrapper from SciPy (look for `*.pyf.src` files in https://github.com/scipy/scipy/tree/main/scipy/linalg)
- add an import to `benchmark/pybench/openblas_wrap/__init__.py`


## Running benchmarks locally

This benchmarking layer is orchestrated from python, therefore you'll need to
have all what it takes to build OpenBLAS from source, plus `python` and

```
$ python -mpip install numpy meson ninja pytest pytest-benchmark
```

The benchmark syntax is consistent with that of `pytest-benchmark` framework. The incantation to run the suite locally is `$ pytest benchmark/pybench/benchmarks/test_blas.py`.

An ASV compatible benchmark suite is planned but currently not implemented.

