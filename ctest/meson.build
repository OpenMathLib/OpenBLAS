testl1_src = ['c_?blat1.f', 'c_?blas1.c']
testl2_src = ['c_?blat2.f', 'c_?blas2.c', 'c_?2chke.c', 'auxiliary.c', 'constant.c', 'c_xerbla.c']
testl3_src = ['c_?blat3.f', 'c_?blas3.c', 'c_?3chke.c', 'auxiliary.c', 'constant.c', 'c_xerbla.c']
testl3_3m_src = ['c_?blat3_3m.f', 'c_?blas3_3m.c', 'c_?3chke_3m.c', 'auxiliary.c', 'constant.c', 'c_xerbla.c']

_test_input_array = {
    'l1': {
        'base': 'x?blat1',
        'has_dat': false,
        'types': ['s', 'd', 'c', 'z'],
        'sources': testl1_src,
    },
    'l2': {
        'base': 'x?cblat2',
        'has_dat': true,
        'types': ['s', 'd', 'c', 'z'],
        'sources': testl2_src,
        'input_file': '?in2',
    },
    'l3': {
        'base': 'x?cblat3',
        'has_dat': true,
        'types': ['s', 'd', 'c', 'z'],
        'sources': testl3_src,
        'input_file': '?in3',
    },
    'l3_3m': {
        'base': 'x?cblat3_3m',
        'has_dat': true,
        'types': ['c', 'z'],
        'sources': testl3_3m_src,
        'input_file': '?in3_3m',
    },
}

_test_runner = executable('test_runner', sources: ['test_runner.c'], install: false)
ctest_inc = _inc + [include_directories('.')]

foreach lvl : ['l1', 'l2', 'l3', 'l3_3m']
    details = _test_input_array[lvl]

    foreach type : details['types']
        op_name = details['base'].replace('?', type)

        mapped_sources = []
        foreach source : details['sources']
            mapped_sources += source.replace('?', type)
        endforeach

        executable(
            op_name,
            sources: mapped_sources + [config_h],
            link_with: [_openblas],
            dependencies: [dependency('threads')],
            include_directories: ctest_inc,
            c_args: ['-DADD_', '-DCBLAS'],
        )

        _args = [f'./@op_name@']
        if details.has_key('input_file')
            _args += [meson.current_source_dir() / details['input_file'].replace('?', type)]
        endif

        test(
            op_name,
            _test_runner,
            args: _args,
            workdir: meson.current_build_dir(),
        ) # TODO: add OPENBLAS_NUM_THREADS=2

    endforeach
endforeach
