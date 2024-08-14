_test_input_array = [
    {'id': 'sblat1', 'has_dat': false},
    {'id': 'sblat2', 'has_dat': true},
    {'id': 'sblat3', 'has_dat': true},
    {'id': 'dblat1', 'has_dat': false},
    {'id': 'dblat2', 'has_dat': true},
    {'id': 'dblat3', 'has_dat': true},
    {'id': 'cblat1', 'has_dat': false},
    {'id': 'cblat2', 'has_dat': true},
    {'id': 'cblat3', 'has_dat': true},
    {'id': 'zblat1', 'has_dat': false},
    {'id': 'zblat2', 'has_dat': true},
    {'id': 'zblat3', 'has_dat': true},
]

if conf_hdat.has('ARCH_X86_64') or conf_hdat.has('ARCH_X86')
    _test_input_array += [
        {'id': 'cblat3_3m', 'has_dat': true}, {'id': 'zblat3_3m', 'has_dat': true},
    ]
endif

_test_runner = executable('test_runner', sources: ['test_runner.c'], install: false)

foreach _test : _test_input_array
    test_id = _test['id']
    executable(
        test_id,
        sources: f'@test_id@.f',
        link_with: [openblas_static],
        dependencies: [dependency('threads')],
    )

    _args = [f'./@test_id@']
    if _test['has_dat']
        _args += [ meson.current_source_dir() / f'@test_id@.dat']
    endif

    test(
        test_id,
        _test_runner,
        args: _args,
        workdir: meson.current_build_dir(),
    )
endforeach
