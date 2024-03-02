# meson setup build --buildtype release
# meson compile -C build
# meson install --prefix=$HOME/.local/lapack

lapack_major_version = 3 # soversion
lapack_minor_version = 12
lapack_patch_version = 0
lapack_version = f'@lapack_major_version@.@lapack_minor_version@.@lapack_patch_version@'

subdir('BLAS/SRC') # Defines blas
subdir('SRC')
subdir('INSTALL')
