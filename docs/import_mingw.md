Microsoft Windows has this thing called "import libraries". You don't need it in MinGW because the `ld` linker from GNU Binutils is smart, but you may still want it for whatever reason.

## Make the `.def`

Import libraries are compiled from a list of what symbols to use, `.def`. This should be already in your `exports` directory: `cd OPENBLAS_TOP_DIR/exports`.

## Making a MinGW import library

MinGW import libraries have the suffix `.a`, same as static libraries. (It's actually more common to do `.dll.a`...)

You need to first prepend `libopenblas.def` with a line `LIBRARY libopenblas.dll`:

    cat <(echo "LIBRARY libopenblas.dll") libopenblas.def > libopenblas.def.1
    mv libopenblas.def.1 libopenblas.def

Now it probably looks like:

    LIBRARY libopenblas.dll
    EXPORTS
	   caxpy=caxpy_  @1
	   caxpy_=caxpy_  @2
           ...

Then, generate the import library: `dlltool -d libopenblas.def -l libopenblas.a`

Again, there is basically **no point** in making an import library for use in MinGW. It actually slows down linking.

## Making a MSVC import library

Unlike MinGW, MSVC absolutely requires an import library. Now the C ABI of MSVC and MinGW are actually identical, so linking is actually okay. (Any incompatibility in the C ABI would be a bug.)

The import libraries of MSVC have the suffix `.lib`. They are generated from a `.def` file using MSVC's `lib.exe`. See [the MSVC instructions](https://github.com/xianyi/OpenBLAS/wiki/How-to-use-OpenBLAS-in-Microsoft-Visual-Studio#generate-import-library-before-0210-version).

## Notes
* Always remember that MinGW is **not the same** as MSYS2 or Cygwin. MSYS2 and Cygwin are full POSIX environments with a lot of magic such as `fork()` and its own `malloc()`. MinGW, which builds on the normal Microsoft C Runtime, has none of that. Be clear about which one you are building for.

