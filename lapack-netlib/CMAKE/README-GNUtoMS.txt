The GNUtoMS CMake module helps LAPACK provide MS-compatible DLLs on
Windows when built with a free GNU Fortran compiler (e.g. MinGW).  If
MS Visual Studio tools are installed when one configures LAPACK to
build with GNU tools the module extends the shared library link rule.
The extended rule creates both a GNU-style .dll.a import library and a
MS-format .lib import library.

LAPACK CMake code installs the import libraries for both formats.
Applications built using CMake can be configured automatically to use
the import libraries matching the target toolchain.
